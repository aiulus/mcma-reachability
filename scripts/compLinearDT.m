% compLinearDT.m
% ------------------------------------------------------------------------------
% This script runs both <enter DDRA ref.> and CORA's conformance pipeline 
% on the same randomly generated dataset. Replace any hardcoded parameters as needed.
%
% Prerequisites (on MATLAB path):
%   • estimateAB_ddra, propagateDDRA    
%   • createDataSet 
%
% Usage:
%   Run this file directly. Adjust “sysType”, “dim”, “initPoints”, “T” to change
%   the system or trajectory lengths.

clear; clc;
%% --- 0.0  Make sure CORA's @polytope wins --------------------------------
% Remove MPT-3 compatibility layer and flush any already-loaded class
%compatRoot = fullfile( ...
%    userpath, 'tbxmanager','toolboxes','mpt','3.2.1','all', ...
%    'mpt3-3_2_1','mpt','modules','compatibility');

%if exist(compatRoot,'dir')
%    rmpath(genpath(compatRoot));                 % drop the whole subtree
%    fprintf('[CORA-fix] removed MPT-3 compatibility folder:\n  %s\n', compatRoot);
%end

% Put CORA’s contSet folder at the very front, just to be safe
%addpath(fullfile('C:\Users\aybuk\Documents\MATLAB\CORA','contSet'),'-begin');

% Flush everything that MATLAB may have cached *before* the path change
%clear classes                                     % unload Polyhedron etc.
%rehash toolboxcache
%% -------------------------------------------------------------------------


%% 0 - Specify system & data parameters
% 'Square': nonlinearARX; 'pedestrian': nonlinearSysDT
%systype = 'mockSysARX'; 
%systype = 'polyNARX';

%systype = 'mockSys';

%% TODO: Defined classes should only contain n=p
%% OR: output-reachset version of DDRA?
systype = 'testSys2';
%dim = 4;

conformance_method = "black";

plot_toggle = struct('ddra', 1, 'cc', 0);

%rng(2);
rand('seed', 1); 

%% 1 - Simulate the system / generate the datasets
% custom_loadDynamics - extends CORA's loadDynamics()
sysparams = struct(); % put x_dim here and pass it to custom_loadDynamics() for scalability analysis
%sysparams.dim = 4;
[sys, params.R0, params.U, params.p_true] = custom_loadDynamics(systype, "rand", sysparams);
dt = sys.dt;

%% (TODO) Consolidate: DDRA models process noise, CC msmt. noise
% uses uncertainty set specifications in loadDynamics, option "standard"

% Create the datasets
% getConfig() - Custom function. Sets hyperparameters such as the number of 
%               distinct trajectories in the testSuite, length of the
%               trajectories, etc. Will later be converted to .mat
%               configuration files. 
cfg = getConfig();
settings = cfg.settings;

%% DEBUG STATEMENT - TODO: Remove when done
settings.n_k = 2;
settings.n_k_train = 2;
settings.n_k_val = 2;
%%

k = settings.n_m_train;
T_k = settings.n_k_train;

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

%% 2 - Run the DDRA pipeline

% aux_CORAtoDDRA - Custom function that converts testSuite objects to data
%                  representation format that the DDRA pipeline expects
[x_all, utraj_all] = narx_CORAtoDDRA(testSuite_train);
%[x_all, utraj_all] = dataTransform_CORA2DDRA(testSuite_train);

% getTrajsDDRA - Custom function. Takes single-trajectory data and creates
%                the time-shifted objects X_-, X_+ etc.
[X_0T, X_1T, U_0T, U_1T, U_full] = shift_trajs(testSuite_train);
%[X_0T, X_1T, U_0T, U_1T, U_full] = shiftTrajs_CORA2DDRA(testSuite_train);

% Initialize data structures for the zonotopes
X0_set = []; U_set = []; 
%W = zonotope(zeros(sys.nrOfOutputs, 1), 0*eye(sys.nrOfOutputs, size(X_1T, 2)));
W = zonotope(zeros(sys.nrOfStates, 1), 0*eye(sys.nrOfStates, size(X_1T, 2)));

% estimateAB_ddra - Custom function. Computes $\mathcal{M}_{AB}$ 
%                   (also annotated as $\mathcal{M}_{\Sigma}$$ according to
%                   Alanwar et.al.
WmatZ = zonotope(0*ones(sys.nrOfOutputs,1),(1e-4)*ones(sys.nrOfOutputs,size(X_1T, 2)));
M_ab = estimateAB_ddra(sys, X_0T, X_1T, U_0T, WmatZ);

totalsteps = settings.n_k_train; % #(identification steps after identification)

% propagateDDRA - Custom function. Uses the previously computed
%                 $\mathcal{M}_{AB}$ to compute the reachable sets. 
[X_model_P2, X_data_P2] = propagateDDRA(U_full, X_0T, X_1T, params.R0, params.U, W, sys, M_ab, totalsteps);


%% 3 - Run the Conformance Checking pipeline
% Pass any relevant parameters to flexBlackBoxConform
sysparams.cfg = cfg;

switch conformance_method
    case "black"
        testSuites_cell = {
            reshapeTestSuiteLinear(testSuite, sys), ...
            reshapeTestSuiteLinear(testSuite_train, sys), ...
            reshapeTestSuiteLinear(testSuite_val, sys)
        };
        [completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, ...
            'testSuites', testSuites_cell, 'sysparams', sysparams);
    case "gray"
        %grayAlg = ["graySim","graySeq","grayLS"]; % pass to
        %flexGrayBoxConform(..., 'grayAlg', grayAlg) - default: graySeq
        [completed, results, R_id, R_val] = flexGrayBoxConform('dynamics', systype, ...
            'testSuites', testSuites, 'sysparams', sysparams);
    case "white"
        [completed, results, R_id, R_val] = flexWhiteBoxConform('dynamics', systype, ...
            'testSuites', testSuites, 'sysparams', sysparams);
end


cc_reachsets = R_id{1};
cc_reachsets = cc_reachsets{1};
ddra_reachsets = X_data_P2;
true_reachsets = X_model_P2;

cc_reachset_norms = [];
for i=1:length(cc_reachsets)
    cc_reachset_norms(i) = norm(cc_reachsets{i}); 
end
ddra_reachset_norms = [];
true_reachset_norms = [];
for i=1:length(ddra_reachsets)
    ddra_reachset_norms(i) = norm(ddra_reachsets{i}); 
    true_reachset_norms(i) = norm(true_reachsets{i});
end


%% --- Plots ---------------------------------------------
fig1 = figure('Name','Zonotope norms over time');
plot(cc_reachset_norms,'-o','LineWidth',1.5);  hold on;
plot(ddra_reachset_norms,'-s','LineWidth',1.5);
plot(true_reachset_norms,'-+','LineWidth',1.5);
grid on;  xlabel('time index k');  ylabel('‖Z_k‖₂');
legend({'CC','DDRA'},'Location','best');
title('Evolution of reach-set norms');
grid on;

outputDir = '../outputs/figures';
if ~exist(outputDir, 'dir')
   mkdir(outputDir)
end
saveas(fig1, fullfile(outputDir, 'reachset_norms_evolution.png'));
fprintf('Saved figure to %s\n', fullfile(outputDir, 'reachset_norms_evolution.png'));


%% TODO: 'visualizeAlinearDT' needs fixing
% Visualize
if plot_toggle.ddra
    projectedDims = {[1 2]};
    numberofplots = length(X_model_P2); %length(X_model_P2)
    visualizeReachsets(params.R0, X_model_P2, X_data_P2, cc_reachsets, projectedDims, numberofplots);
    %visualizeDDRA(params.R0, X_model_P2, X_data_P2, ...
    %          projectedDims, axx, numberofplots);
end