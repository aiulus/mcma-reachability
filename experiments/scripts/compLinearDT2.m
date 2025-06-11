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
%systype = 'chain_of_integrators';
systype = 'Square';
dim = 4;
dt = 0.05;

initpoints = 5; % #(distinct trajectories to simulate)
T = 120; % length of each trajectory

plot_toggle = struct('ddra', 0, 'cc', 0);

rng(2);

%% 1 - Simulate the system / generate the datasets
[sys, params_true.R0, params_true.U, p_true] = custom_loadDynamics(systype, "rand", dim);
cfg = getConfig();
settings = cfg.settings;
n_k_total = settings.n_k + settings.n_k_train + settings.n_k_val;
n_m_total = settings.n_m + settings.n_m_train + settings.n_m_val;
n_s_total = settings.n_s + settings.n_s_train + settings.n_s_val;
testSuite = createTestSuite(sys, params_true, n_k_total, n_m_total, n_s_total, cfg.options_testS);

% Initialize data structures for the zonotopes
X0_set = []; U_set = []; W = []; WmatZ = [];

[X0, U, W, Wmatzono] = initialSetupDDRA(sys, initpoints, T, ...
                                       ones(dim,1), 0.1*eye(dim), ...   % X0_center & X0_spread
                                       1*ones(sys.dims.m,1), 0.25*eye(sys.dims.m), ...   % U_center & U_spread
                                       0*ones(dim,1), 0.005*eye(dim));  % W_center & W_spread

% Simulate the system dynamics with random samples from X0 and W
%[x_all, utraj_all] = getDataDDRA(sys, initpoints, T, X0, U, W);

% Create the datasets
%[data, testSuites] = createDataSet(sys, initpoints, T, x_all, utraj_all);

%% 2 - Run the DDRA pipeline
[x_all, utraj_all] = convertDynamicsData(testSuite, sys);
[U_full, X_0T, X_1T] = getTrajsDDRA(sys, initpoints, T, x_all, utraj_all, false);
M_ab = estimateAB_ddra(sys.discrete, X_0T, X_1T, U_full, Wmatzono);
totalsteps = 10; % #(identification steps after identification)
[X_model_P2, X_data_P2] = propagateDDRA(X0, U, W, sys.discrete, M_ab, totalsteps);

% Visualize
if plot_toggle.ddra
    projectedDims = {[1 2],[3 4],[4 5]};
    axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
    numberofplots = length(X_model_P2); %length(X_model_P2)
    visualizeAlinearDT(X0, X_model_P2, X_data_P2, projectedDims, axx, numberofplots);
end

%% 3 - Run the Conformance Checking pipeline
%[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'testSuites', testSuites, 'sysparams', dim);
%% Data-passing temporarily disabled
[completed, results, R_id, R_val] = flexBlackBoxConform('dynamics', systype, 'sysparams', dim);