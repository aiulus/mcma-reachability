% compLinearDT.m
% ------------------------------------------------------------------------------
% This script runs both <enter DDRA ref.> and CORA's black-box ARX
% conformance pipeline on the same randomly generated dataset. Replace any
% hardcoded parameters as needed.
%
% Prerequisites (on MATLAB path):
%   • systemsDDRA, initialSetupDDRA, getDataDDRA, getTrajsDDRA,
%     estimateAB_ddra, propagateDDRA      (Paper 2)
%   • example_nonlinearARX_conform_03_black_common,
%     convertCommonToP1                    (Paper 1)
%   • createDataSet (aux_DDRA + aux_CC)    (this file)
%
% Usage:
%   Run this file directly. Adjust “sysType”, “dim”, “initPoints”, “T” to change
%   the system or trajectory lengths.

clear; clc;

%% 0 - Specify system & data parameters
% 'Square': nonlinearARX; 'pedestrian': nonlinearSysDT
systype = 'lipschitzNARX'; 
%systype = 'polyNARX';
dim = 4;
dt = 0.1;

plot_toggle = struct('ddra', 0, 'cc', 0);

rng(2);

%% 1 - Simulate the system / generate the datasets
% custom_loadDynamics - extends CORA's loadDynamics()
sysparams.dim = dim;
[sys, params.R0, params.U, params.p_true] = custom_loadDynamics(systype, "rand", sysparams);


%% (TODO) Consolidate: DDRA models process noise, CC msmt. noise
% uses uncertainty set specifications in loadDynamics, option "standard"

% initialSetupDDRA - just sets the uncertainty sets
%[~, ~, W, Wmatzono] = initialSetupDDRA(sys, initpoints, T, ...
%                                       0, 0, ...   % X0_center & X0_spread
%                                       0.1, 0.2, ...   % U_center & U_spread
%                                       -0.05, 0.1);  % W_center & W_spread

% Create the datasets
% getConfig() - Custom function. Sets hyperparameters such as the number of 
%               distinct trajectories in the testSuite, length of the
%               trajectories, etc. Will later be converted to .mat
%               configuration files. 
cfg = getConfig();
settings = cfg.settings;
%settings.n_s = 1;
%settings.n_s_train = 1;
%settings.n_s_val = 1;
%params = struct('R0', R0, 'U', U);
%initpoints = 5; % #(distinct trajectories to simulate)
%T = 120; % length of each trajectory
k = settings.n_m_train;
T_k =settings.n_k_train;

% createTestSuite - CORA fuction
testSuite = createTestSuite(sys, params, settings.n_k, settings.n_m, settings.n_s, cfg.options_testS);
testSuite_train = createTestSuite(sys, params, settings.n_k_train, settings.n_m_train, settings.n_s_train);
testSuite_val = createTestSuite(sys, params, settings.n_k_val, settings.n_m_val, settings.n_s_val);

testSuites = cell(3,1);
testSuites{1} = testSuite;
testSuites{2} = testSuite_train;
testSuites{3} = testSuite_val;

% union_testSuites - Custom function. Builds the union of multiple
%                    testSuite objects. 
%complete_testSuite = union_testSuites(testSuite, testSuite_train, testSuite_val);

%% ACHTUNG!!-- aux_CORAtoDDRA currently propagates X0 with (A, B, C, D)
%%          -- only use systems with n=p in the future!

% aux_CORAtoDDRA - Custom function that converts testSuite objects to data
%                  representation format that the DDRA pipeline expects
%[x_all, utraj_all] = narx_CORAtoDDRA(complete_testSuite, sys);
[x_all, utraj_all] = lsdt_CORAtoDDRA(testSuite_train);

%% 2 - Run the DDRA pipeline
% getTrajsDDRA - Custom function. Takes single-trajectory data and creates
%                the time-shifted objects X_-, X_+ etc.
[X_0T, X_1T, U_plus] = shift_trajs(T_k, x_all, utraj_all);

% Initialize data structures for the zonotopes
X0_set = []; U_set = []; 
W = zonotope(zeros(sys.nrOfStates, 1), 0*eye(sys.nrOfStates, size(X_1T, 2)));
WmatZ = zonotope(zeros(sys.nrOfStates, 1), 0.01*eye(sys.nrOfStates, size(X_1T, 2)));

% estimateAB_ddra - Custom function. Computes $\mathcal{M}_{AB}$ 
%                   (also annotated as $\mathcal{M}_{\Sigma}$$ according to
%                   Alanwar et.al.
M_ab = estimateAB_ddra(sys, X_0T, X_1T, U_plus, WmatZ);

totalsteps = 10; % #(identification steps after identification)

% propagateDDRA - Custom function. Uses the previously computed
%                 $\mathcal{M}_{AB}$ to compute the reachable sets. 
[X_model_P2, X_data_P2] = propagateDDRA(params.R0, params.U, W, sys, M_ab, totalsteps);

% Visualize
if plot_toggle.ddra
    projectedDims = {[1 2],[3 4],[4 5]};
    axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
    numberofplots = length(X_model_P2); %length(X_model_P2)
    visualizeAlinearDT(X0, X_model_P2, X_data_P2, projectedDims, axx, numberofplots);
end

%% 3 - Run the Conformance Checking pipeline
% Pass any relevant parameters to flexBlackBoxConform
sysparams.cfg = cfg;
[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'testSuites', testSuites, 'sysparams', sysparams);

% Temporarily disabling dataset passing
%[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'sysparams', sysparams);