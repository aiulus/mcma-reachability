function completed = example_chainOfIntegrators
% example_chainOfIntegrators -

% ------------------------------ BEGIN CODE -------------------------------

%% Initial Setup

%dyn = "chain_of_integrators";
dyn = "pedestrian";
dt = 0.5;

n_n = 2; % Order of the integrator chain

n_k = 8; % number of time steps for identification
n_k_val = 20; % number of test cases for identification
n_m = 20; % number of test cases for identification
n_m_val = 100; % number of test cases for validation

n_s = 1; % number of samples per test case

constraints = {'gen','half',}; % type of containment constraints

rng(2)

% Reachability sSettings
options_reach.zonotopeOrder = 100;
options.reach = options_reach;

% Conformance Settings
options.cs.robostnessMargin = 1e-9;
options.cs.cost = 'interval';
options.cs.verbose = false;
options_testS.p_extr = 0.2; % Probability of extreme points
options_testS.inputCurve = 'rand'; % Input trajectory shape

% Evaluation settings
check_contain = false;
plot_settings.dims = [1, 2];
plot_settings.name = sprintf('Conformance Synthesis: %s', dyn, n_n);

% Parameters and System Dynamics
if dyn == "chain_of_integrators"
    [sys, params_true.R0, params_true.U] = extended_loadDynamics(dyn, dt, n_n);
else 
    [sys, params_true.R0, params_true.U] = loadDynamics(dyn); % TODO: Add list of options for this to the summary
end
params_true.tFinal = dt * (n_k - 1);

% Simulation
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, options_testS);

% Initial estimates of the disturbance sets
c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0, eye(sys.dims.n)]);
params_id_init.U = zonotope([c_U, eye(sys.dims.m)]);

%% Conformance identification
num_id = length(constraints);

% Struct for saving identification results for each system
%% (Q) for each constraint?
configs = cell(num_id + 1, 1);
configs{1}.sys = sys;
configs{1}.params = params_true;
configs{1}.options = options_reach;
configs{1}.name = 'true';

for i_id = 1:num_id
    options.cs.constraints = constraints{i_id};
    %% (Q) What's with the 3*n_n here?
    fprintf(['Identification with %s-constraints, n_m = %d, n_k = %d, ' ...
        'n_x = %d\n'], options.cs.constraints, n_m, n_k, 3*n_n);
    tic;
    [configs{i_id + 1}.params, ~] = conform(sys, params_id_init, options);
    configs{i_id + 1}.sys = sys;
    configs{i_id + 1}.options = options_reach;
    configs{i_id + 1}.name = options.cs.constraints;
    Ts = toc;
    fprintf('Identification time: %.4f\n', Ts);
end

%% Validation and Visualization

% Create Validation Data
if n_m_val ~= 0
    params_true.tFinal = dt * (n_k_val - 1);
    testSuite_val = createTestSuite(sys, params_true, n_k_val, ...
        n_m_val, n_s, options);
end

% Combine validation and training test cases
testSuite{1} = combineTestCases(params_true.testSuite{1}, testSuite_val{1});
plot_settings.s_val = size(params_true.testSuite{1}.y, 3) + 1;

validateReach(testSuite{1}, configs, check_contain, plot_settings);

completed = true;

end