function out = func_ddra_linearDT(systype, dim, scale_noise, plot_toggle)
% runDDRA_linear - Execute matrix zonotope-based reachability analysis
% on discrete-time linear systems.
%
% Inputs:
%   system_type  - (char) system type name (e.g., 'chain_integrators')
%   dim          - (int) system dimensionality (only for scalable types)
%   scale_noise  - scaling for process noise
%   plot_toggle  - (logical) plot trajectory and reachable set projections

%------------- BEGIN CODE --------------
rng(1);

%% Define the system
dt = 0.05;
sys = systemsDDRA(systype, dt, dim);
sys_d = sys.discrete;
% sys_c = sys.cont;
n = sys.dims.n;

%% Initial setup
initpoints = 1; % #trajectories
steps = 120; % #time steps

X0_center = 1; X0_spread = 0.1;
U_center = 10; U_spread = 0.25;
W_center = 0; W_spread = 0.005 * scale_noise;


[X0, U, W, Wmatzono] = initialSetupDDRA(sys, initpoints, steps, ...
    X0_center, X0_spread, U_center, U_spread, W_center, W_spread);


%% Simulate & get the data
[x, utraj] = getDataDDRA(sys, initpoints, steps, X0, U, W);

%% Get the trajectories
[U_full, X_0T, X_1T] = getTrajsDDRA(sys, initpoints, steps, x, utraj, plot_toggle);
    

%% Estimate M_ab
M_ab = estimateAB_ddra(sys_d, X_0T, X_1T, U_full, Wmatzono);


%% compute next step sets from model / data
totalsteps = 5; % #propagation steps
[X_model, X_data] = propagateDDRA(X0, U, W, sys_d, M_ab, totalsteps);

%% Output
out = struct('X_data', {X_data}, 'X_model', {X_model});

%% Visualize
if plot_toggle
    projectedDims = {[1 2],[3 4],[4 5]};
    axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
    numberofplots = 5; %length(X_model)
    visualizeAlinearDT(X0, X_model, X_data, projectedDims, axx, numberofplots);
end
%------------- END OF CODE --------------



