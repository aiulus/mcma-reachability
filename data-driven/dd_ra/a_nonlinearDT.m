% t_nonlinearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = f(x(k),u(k)) + w(k)  
% The approach is based on [1]. The nonlinear system is found in [2]
% 
% 
%
% Syntax:  
%    t_nonlinearDT
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% Example: 
%
% References:
%    [1] Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson 
%       "Data Driven Reachability Analysis Using Matrix Zonotopes"
%    [2] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006
%    
% 
% Author:       Amr Alanwar
% Written:      29-October-2020
%
% Adapted by:   Aybüke Ulusarslan
% Updated:      25-April-2025
% Last update:  
% Last revision:---


%------------- BEGIN CODE --------------
rng(1);
clearvars;


%% Initial setup
% System initialization
systype = 'cstr';
dt = 0.015;
sys = coraSystemsWrapper(systype, dt);

% Centers and spread factors for the zonotopes
u0 = [0.01;0.01]; Gu0 = diag([0.1;0.2]);
y0 = [-1.9;-20]; Gy0 = diag([0.005;0.3]);
w0 = 0; Gw = 1e-4;

% Catchall dictionary
lookup = struct( ...
    'sys', sys, ...
    'dt', dt, ...
    'NN', 5, ...
    'initpoints', 30, ...
    'steps', 20, ...
    'U', zonotope(u0, Gu0), ... % Input set
    'R0', zonotope(y0, Gy0), ... % Initial reachable set
    'W', zonotope(w0*ones(sys.dims.n,1),Gw*ones(sys.dims.n,1)), ... % Measurement noise bound
    'n', sys.dims.n, ... % Input dimension
    'fun', sys.fun, ... % System dynamics
    'zonotopeOrder', 100, ... % Reachability settings 
    'tensorOrder', 2, ... % Reachability settings 
    'errorOrder', 5 ... % Reachability settings 
    );

lookup.tFinal = lookup.dt*lookup.NN;
lookup.totalsamples = lookup.steps*lookup.initpoints; % total #samples


Wmatzono= getGW(lookup);

%% Allocate trajectories
[x_meas_vec_0, x_meas_vec_1, x_free_vec_0, x_free_vec_1, U_full, X_0T, X_1T] = getTrajsNonlinDDRA(lookup);
lookup.U_full = U_full;
lookup.X_0T = X_0T;
lookup.X_1T = X_1T;

%% ?
stepsLip = 1;
initpointsLip = 50;
[gamma,L] = compLipConst(lookup.fun, lookup.U, lookup.R0, stepsLip, initpointsLip, lookup.n);
eps(1) = L(1) .* gamma(1)/2;
eps(2) = L(2) .* gamma(2)/2;
lookup.Zeps = zonotope([zeros(2,1),diag(eps)]);
lookup.ZepsFlag = 1;

%% Define the discrete system 
sysDisc = nonlinearDT('stirredTankReactor', @(x,u) cstrDiscr(x,u,dt), dt, sys.dims.n, sys.dims.m);

params.R0 = lookup.R0;
params.U = lookup.U;

%% Reachability Analysis 
% compute model based reachability (R) and data driven one (R_data)
tic
[R ,R_data]= reach_DT(sysDisc, params, lookup);
tComp = toc;
disp("Computation time: " + tComp);

if lookup.ZepsFlag
    for i = 1:NN
        R_data.timePoint.set{i} = R_data.timePoint.set{i} + lookup.Zeps;
    end
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot initial set
handleX0=plot(R0,[1,2],'k-','LineWidth',2);

% plot model based reachable set
handleModel=plot(R,[1 2],'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');

% plot data driven reachable set
handleData=plot(R_data,[1 2],'r','Filled',false);


% formatting
xlabel('x_1');
ylabel('x_2');

% skip warning for extra legend entries
warOrig = warning; warning('off','all');
legend([handleX0,handleModel,handleData],...
    'Initial set $\mathcal{X}_0$','Set from model $\mathcal{R}_k$','Set from data $\mathcal{R}_k^{\prime}$','Location','northwest','Interpreter','latex');
warning(warOrig);
ax = gca;
ax.FontSize = 16;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% example completed
completed = 1;

%------------- END OF CODE --------------