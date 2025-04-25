% t_linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) + w(k)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis Using Matrix Zonotopes"
%
%
%
% Author:       Amr Alanwar
% Written:      28-October-2020
%
% Adapted by:   Aybüke Ulusarslan
% Edited:       25-April-2025
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------
rng(1);
clearvars;

%% Define the system
dt = 0.05;
sys = systemsDDRA('example0', dt);
sys_d = sys.discrete;
sys_c = sys.cont;
n = sys.dims.n;

%% Initial setup
initpoints = 1; % #trajectories
steps = 120; % #time steps

X0_center = ones(n,1); X0_spread = 0.1;
U_center = 10; U_spread = 0.25;
W_center = 0; W_spread = 0.005;


[X0, U, W, Wmatzono] = initialSetupDDRA(sys, initpoints, steps, ...
    X0_center, X0_spread, U_center, U_spread, W_center, W_spread);


%% Simulate & get the data
[x, utraj] = getDataDDRA(sys, initpoints, steps, X0, U, W);

%% Get the trajectories
plot_toggle = 1;
[u_mean_vec_0, x_meas_vec_0, x_meas_vec_1, U_full, X_0T, X_1T] = ...
    getTrajsDDRA(sys, initpoints, steps, x, utraj, plot_toggle);

%% Estimate M_ab
[M_ab, intAB1] = estimateAB_ddra(sys_d, X_0T, X_1T, U_full, Wmatzono);


%% compute next step sets from model / data
totalsteps = 5; % #propagation steps
[X_model, X_data] = propagateDDRA(X0, U, W, sys_d, M_ab, totalsteps);


%% Visualize
projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
numberofplots = 5; %length(X_model)
visualizeAlinearDT(X0, X_model, X_data, projectedDims, axx, numberofplots);

%------------- END OF CODE --------------

