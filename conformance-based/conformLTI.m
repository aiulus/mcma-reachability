function [runtimes, Rkx] = conformLTI(systype, x_dim, unc_scaling, optionals)
% ADAPTED FROM:
%   - CORA/examples/contDynamics/linearSysDT/example_linearSysDT_conform_04_LTV.m,
%   - CORA/examples/contDynamics/nonlinearARX/example_nonlinearARX_conform_03_black
%% User Specifications ----------------------------------------------------

%% TODO:
% 1- unc_scaling should be used to play around with the uncertainty sets,
%    primarily defined by loadDynamics.m (CORA)
% 2- Rkx should store the reachsets

dyn = systype; % dynamics (choose from "platoon", "pedestrian")

% Set Random Number Stream
rng(2)

% Get default settings & Override with 'optionals' if appicable
settings = getStandardSettings();

% cost_norm = "interval"; % Reachable set size evaluation method
constraints = "half"; % Halfspace constraints for containment
methodsGray = ["blackGP", "blackCGP"];

if nargin >= 4 && ~isempty(optionals)
    fields = fieldnames(optionals);
    for f = 1:numel(fields)
        settings.(fields{f}) = optionals.(fields{f});
    end
end

% Evaluation settings
check_contain = false;
plot_settings.dims = [1 2]; % TODO: Project to couples of axes for x_dim > 2 & change to plot_settings.dims = 1:x_dim
plot_settings.name = sprintf("Conformance Synthesis: %s", dyn);

% Parameter extraction
n_k = settings.n_k;
n_k_val = settings.n_k_val;
n_m = settings.n_m;
n_m_val = settings.n_m_val;
n_s = settings.n_s;
constraints = settings.constraints;
options_testS = settings.options_testS;
options_reach = settings.options_reach;

%% TODO: Remove hard-coding
options.approx.p = 3; % Regression order of the ARX model
options.approx.verbose = 1;

% Parameters and System Dynamics
if dyn == "platoon"
    n_n = x_dim;
    [sys, params_true.R0, params_true.U] = aux_load_platoon(n_n,...
        max(n_k,n_k_val), unc_scaling);
else
    [sys, params_true.R0, params_true.U] = custom_loadDynamics(dyn, unc_scaling);
end
params_true.tFinal = sys.dt * n_k - sys.dt;

% Initial Estimates of the Disturbance Sets
c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfStates)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs)]);

% Simulation
params_id_init.testSuite_train = createTestSuite(sys, params_true, n_k, n_m, n_s, options_testS);
%params_id_init.testSuite_train = createTestSuite(sys, params_true, 2, 13, 19, options_testS);
params_id_init.testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s, options_testS);

%% Conformance Identification ---------------------------------------------
num_id = length(constraints);
%name_id = cell(num_id,2);

% Struct for saving the identification results for each system
configs = cell(num_id+1,1);
configs{1}.sys = sys;
configs{1}.params = params_true;
configs{1}.options = options_reach;
configs{1}.name = "true";

for i_id = 1:num_id
    % run the identification
    options.cs.constraints = constraints{i_id};
    fprintf("Identification with %s-constraints, n_m=%d, " + ...
        "n_k=%d, n_x=%d\n",options.cs.constraints, n_m, n_k, 3*n_n);
    tic;
    [configs{i_id+1}.params, results] = custom_conform_black(params_id_init, options, "blackCGP");
    configs{i_id+1}.sys = sys;
    configs{i_id+1}.options = options_reach;
    configs{i_id+1}.name = options.cs.constraints;
    Ts=toc;
    fprintf("Identification time: %.4f\n", Ts);
end

%% Validation and Visualization -------------------------------------------

% Create Validation Data
if n_m_val ~= 0
    params_true.tFinal = sys.dt * n_k_val - sys.dt;
    testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
        n_s, options);
end
% combine validation and trainings test cases
testSuite{1} = combineTestCases(params_true.testSuite{1}, testSuite_val{1});
plot_settings.s_val = size(params_true.testSuite{1}.y,3) + 1; 
    % (setting s_val leads to different color for validation test cases)

% run validation and plotting
validateReach(testSuite{1}, configs, check_contain, plot_settings);

% example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys, R0, U] = aux_load_platoon(N_v,N_k, unc_scaling)
% load the tank dynamics with the specified dimension--
% called with n_n: #vehicles
%             max(...
%                 n_k: #time steps for identification, 
%             n_k_val: #time test cases for idetification
%                 )

dt = 0.5; % Time step for discretization
N_u = N_v; % #vehicles
N_n = N_v*3; % 3*#vehicles
sys = platoonN(dt,N_v,N_k);

c_R0 = randn(N_n,1); 
for i=0:N_v-1
    c_R0(i*N_v + 2) = 3*abs(c_R0(i*N_v + 2));
end
alpha_R0 = unc_scaling.R0 * 2 * rand(N_n,1);
c_U = randn(N_u,1);
alpha_U = unc_scaling.U * rand(N_u,1);
R0 = zonotope([c_R0,diag(alpha_R0)]);
U = zonotope([c_U,diag(alpha_U)]);

end

function settings = getStandardSettings 
    % Conformance identification
    settings.n_m = 2; % #different input trajectories
    settings.n_s = 50; % #sample trajectories per input trajectory
    settings.n_k = 4; % length (#time steps) of identification trajectories

    % Training and validation data
    settings.n_m_train = 100;
    settings.n_s_train = 10;
    settings.n_k_train = 4;
    settings.n_m_val = 5;
    settings.n_k_val = 4;       
   
    % Reachability Settings
    settings.options_reach.zonotopeOrder = 100;
    settings.options_reach.tensorOrder = 2;
    settings.options_reach.errorOrder = 1;
    settings.options_reach.tensorOrderOutput = 2;
    settings.options_reach.verbose = false;

    % testSuite settings
    settings.options_testS.p_extr = 0.3;
    
    % Conformance Settings
    settings.options = settings.options_reach;
    settings.options.cs.robustnessMargin = 1e-9;
    settings.options.cs.cost = 'interval';
    settings.options.cs.verbose = false;
    
    % Black-box approximation options
    settings.options.approx.gp_parallel = true;
    settings.options.approx.gp_pop_size = 50;
    settings.options.approx.gp_num_gen = 30;
    settings.options.approx.gp_func_names = {'times','plus', 'square'};
    settings.options.approx.gp_max_genes = 2;
    settings.options.approx.gp_max_depth = 2;
    settings.options.approx.gp_parallel = false;
    settings.options.approx.cgp_num_gen = 5;
    settings.options.approx.cgp_pop_size_base = 5;
    settings.options.approx.save_res = false;
    settings.options.approx.p = sys.n_p;
end

% ------------------------------ END OF CODE ------------------------------

