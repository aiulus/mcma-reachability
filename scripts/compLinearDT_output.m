%% 0 - Specify system & data parameters

%% TODO: Defined classes should only contain n=p
%% OR: output-reachset version of DDRA?
systype = 'testSys2';
%systype = 'nonlinearARX';
%dim = 4;

conformance_method = "gray";

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

T_k = settings.n_k_train;
steps = T_k;
initpoints = settings.n_m_train;
n_s = settings.n_s_train;

% createTestSuite - CORA fuction
suite_test = createTestSuite(sys, params, settings.n_k, settings.n_m, settings.n_s, cfg.options_testS);
suite_train = createTestSuite(sys, params, settings.n_k_train, settings.n_m_train, settings.n_s_train);
suite_val = createTestSuite(sys, params, settings.n_k_val, settings.n_m_val, settings.n_s_val);

testSuites = cell(3,1);
testSuites{1} = suite_test;
testSuites{2} = suite_train;
testSuites{3} = suite_val;

% union_testSuites - Custom function. Builds the union of multiple
%                    testSuite objects. 
complete_testSuite = union_testSuites(suite_test, suite_train, suite_val);

%% ACHTUNG!!-- aux_CORAtoDDRA currently propagates X0 with (A, B, C, D)
%%          -- only use systems with n=p in the future!

%% 2 - Run the DDRA pipeline

% aux_CORAtoDDRA - Custom function that converts testSuite objects to data
%                  representation format that the DDRA pipeline expects
%[x_all, utraj_all] = narx_CORAtoDDRA(testSuite_train);
[x_all, utraj_all] = narx_CORAtoDDRA(suite_train);
%[x_all, utraj_all] = dataTransform_CORA2DDRA(testSuite_train);

% getTrajsDDRA - Custom function. Takes single-trajectory data and creates
%                the time-shifted objects X_-, X_+ etc.
[Y_0T, Y_1T, U_0T, U_1T, U_full] = shift_trajs(suite_train);
%[X_0T, X_1T, U_0T, U_1T, U_full] = shiftTrajs_CORA2DDRA(testSuite_train);


% estimateAB_ddra - Custom function. Computes $\mathcal{M}_{AB}$ 
%                   (also annotated as $\mathcal{M}_{\Sigma}$$ according to
%                   Alanwar et.al.
x0GenOrder = 1; 
uGenOrder = 1;
[X0, U, W, M_w, V, M_v, AV, M_Av] = initialSetupDDRA_output( ...
    sys, initpoints, steps, n_s, ...
    zeros(sys.nrOfStates, 1), 0.01, ... X0
    ones(sys.nrOfInputs, 1), 0.001, ... % U
    zeros(sys.nrOfStates, 1), 0.01, ... % Process noise
    zeros(sys.nrOfStates, 1), 0.02, ... % Measurement noise
    x0GenOrder, uGenOrder);
M_sigma = estimateAB_output(sys, Y_0T, Y_1T, U_0T, M_v, M_w, M_Av);


totalsteps = settings.n_k_train; % #(identification steps after identification)

% propagateDDRA - Custom function. Uses the previously computed
%                 $\mathcal{M}_{AB}$ to compute the reachable sets. 
[R_tr_ddra, R_pred_ddra] = propDDRA_output(sys, utraj, X0, W, V, AV, M_sigma, timesteps);

%% 3 - Run the Conformance Checking pipeline
% Pass any relevant parameters to flexBlackBoxConform
sysparams.cfg = cfg;

switch conformance_method
    case "black"
        testSuites_cell = {
            reshapeTestSuiteLinear(suite_test, sys), ...
            reshapeTestSuiteLinear(suite_train, sys), ...
            reshapeTestSuiteLinear(suite_val, sys)
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
ddra_reachsets = R_pred_ddra;
true_reachsets = R_tr_ddra;

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
    numberofplots = length(R_tr_ddra); %length(X_model_P2)
    visualizeReachsets(params.R0, R_tr_ddra, R_pred_ddra, cc_reachsets, projectedDims, numberofplots);
    %visualizeDDRA(params.R0, X_model_P2, X_data_P2, ...
    %          projectedDims, axx, numberofplots);
end