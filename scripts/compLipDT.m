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
systype = 'lipschitzSysDT'; 
%systype = 'polyNARX';
dim = 4;

conformance_method = 'gray';

plot_toggle = struct('ddra', 0, 'cc', 0);

rng(2);

%% 1 - Simulate the system / generate the datasets
% custom_loadDynamics - extends CORA's loadDynamics()
sysparams.dim = dim;
[sys, params.R0, params.U, params.p_true] = custom_loadDynamics(systype, "rand", sysparams);
dt = sys.dt;

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
[x_all, utraj_all] = dataTransform_CORA2DDRA(testSuite_train);

%% 2 - Run the DDRA pipeline
% getTrajsDDRA - Custom function. Takes single-trajectory data and creates
%                the time-shifted objects X_-, X_+ etc.
[X_0T, X_1T, U_0T, U_1T] = shiftTrajs_CORA2DDRA(testSuite_train);

% Initialize data structures for the zonotopes
X0_set = []; U_set = []; 
W = zonotope(zeros(sys.nrOfOutputs, 1), 0*eye(sys.nrOfOutputs, size(X_1T, 2)));
%WmatZ = zonotope(zeros(sys.nrOfOutputs, 1), 0.01*eye(sys.nrOfOutputs, size(X_1T, 2)));
%WmatZ = matZonotope(zeros(sys.nrOfOutputs, (size(X_1T, 2))), {zeros(sys.nrOfOutputs, size(X_1T, 2))});
WmatZ = zonotope(0*ones(sys.nrOfOutputs,1),(1e-4)*ones(sys.nrOfOutputs,1));

% estimateAB_ddra - Custom function. Computes $\mathcal{M}_{AB}$ 
%                   (also annotated as $\mathcal{M}_{\Sigma}$$ according to
%                   Alanwar et.al.
%M_ab = estimateAB_ddra(sys, X_0T, X_1T, U_0T, WmatZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lookup = struct( ...
    'fun', sys.mFile, ...
    'U', params.U, ...
    'R0', params.R0, ...
    'n', sys.nrOfOutputs, ...
    'tensorOrder', cfg.options_reach.tensorOrder ...
    );
stepsLip = 1;
initpointsLip = 50;
[gamma,L] = compLipConst(lookup.fun, lookup.U, lookup.R0, stepsLip, initpointsLip, lookup.n);
eps(1) = L(1) .* gamma(1)/2;
eps(2) = L(2) .* gamma(2)/2;
lookup.Zeps = zonotope([zeros(2,1),diag(eps)]);
lookup.ZepsFlag = 1;

%% Define the discrete system 
options = struct( ...
                'R0', params.R0, ...
                'U', params.U, ...
                'tStart', 0, ...
                'tFinal', 10, ...
                'zonotopeOrder', cfg.options_reach.zonotopeOrder, ...
                'tensorOrder', cfg.options_reach.tensorOrder ...
                );

[X_model_P2 ,X_data_P2]= reach_DT(sys, params, options);
tComp = toc;

if lookup.ZepsFlag
    for i = 1:NN
        X_data_P2.timePoint.set{i} = X_data_P2.timePoint.set{i} + lookup.Zeps;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%totalsteps = 10; % #(identification steps after identification)

% propagateDDRA - Custom function. Uses the previously computed
%                 $\mathcal{M}_{AB}$ to compute the reachable sets. 
%[X_model_P2, X_data_P2] = propagateDDRA(params.R0, params.U, W, sys, M_ab, totalsteps);

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
%[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'testSuites', testSuites, 'sysparams', sysparams);

% Temporarily disabling dataset passing
%[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'sysparams', sysparams);

switch conformance_method
    case "black"
        [completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, ...
            'testSuites', testSuites, 'sysparams', sysparams);
    case "gray"
        %grayAlg = ["graySim","graySeq","grayLS"]; % pass to
        %flexGrayBoxConform(..., 'grayAlg', grayAlg) - default: graySeq
        [completed, results, R_id, R_val] = flexGrayBoxConform('dynamics', systype, ...
            'testSuites', testSuites, 'sysparams', sysparams);
    case "white"
        [completed, results, R_id, R_val] = flexWhiteBoxConform('dynamics', systype, ...
            'testSuites', testSuites, 'sysparams', sysparams);
end