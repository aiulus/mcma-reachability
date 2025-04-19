% ZPC: run ZPC with model and without the model
% i.e., run 
% 1-Robust data driven Predictive control scheme (ZPC) 
% 2- Same ZPC while knowing the model (RMPC-zono)
%
% Inputs:
%    none
%
% Outputs:
%    saved workspace
%
% Example: 
%
% See also: ---

% Original Author:       Amr Alanwar, Yvonne Stürz 
% Adapted by:            Aybüke Ulusarslan, Yongkuan Zhang
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%% INITIAL SETUP
rand('seed',4500);

clear all
close all

%% SYSTEM SETUP
systype = 'example0'; % system type
dt = 0.05; % sampling time
sys = systemsZonoDDSF(systype, dt);
sys_c = sys.cont;
sys_d = sys.discrete;
n = sys.dims.n;

%% SIMULATION SETUP
initpoints =100; % number of trajectories
steps =5; % number of steps for each trajectory
totalsamples = initpoints*steps; %Total number of samples

%% RETRIEVE BOUNDARY CONDITIONS
sys = setupBoundaryConditions(sys);
intc = sys.bcs.intc;
y0 = sys.bcs.y0;
X0 = sys.bcs.X0;
U = sys.bcs.U;

%% SET UP NOISE TERMS
w_spread = 0.01; % controlls process noise levels
v_spread = 0.002; % controlls measurement noise levels
noise_catchall = setupInitialNoises(sys, w_spread, v_spread, totalsamples);

%%% RETRIEVE PROCESS NOISE DETAILS
W = noise_catchall.process.W;
GW = noise_catchall.process.GW;
GmatW = noise_catchall.process.GmatW;
WmatZono = noise_catchall.process.WmatZono;

%%% RETRIEVE MEASUREMENT NOISE DETAILS
V = noise_catchall.msmt.V;
GV = noise_catchall.msmt.GV;
GmatV = noise_catchall.msmt.GmatV;
VmatZono = noise_catchall.msmt.VmatZono;

%%% RETRIEVE PROJECTED MEASUREMENT NOISE DETAILS
AV = noise_catchall.Av.AV;
VAmatZono = noise_catchall.Av.VAmatZono;

%% generate data
[u_traj, x_v, x] = genData(sys, X0, U, W, V, initpoints, totalsamples);

%% prepeare Y_+ Y_-
[U_data, Y_0T, Y_1T] = data2vec(u_traj, x_v, x, initpoints, n, steps);

% plot simulated trajectory
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;

% prepare M_Sigma which is a set of [A B]
M_Sigma = estimateAB(Y_0T, U_data, Y_1T, VmatZono, WmatZono, VAmatZono, sys_d);


%% Compute ZPC problem
%Horizon N for ZPC
N = 2;
% define output cost matrix
%output_cost_coeff = 1e3;
%Qy = output_cost_coeff * eye(sys.dims.n); 
% control cost matrix
input_cost_coeff = 0.001;
Q = input_cost_coeff * eye(sys.dims.m);

% ZPC number of time steps
maxsteps = 80;
% time step for plotting 
timestep_plot = 10;

[uPred, uPred_model, y_t, y_t_model, execTimeZPC, execTimeRMPC] =  runZonoDDSF( ...
    sys_d, y0, intc, r_u, r_y, U, N, maxsteps, timestep_plot, Q, W, V, AV, M_Sigma);          

Cost_model=0;
for i=1:maxsteps
    Cost_model_vec(i) = (uPred_model(:,i) - u_l(:, i))'*Q*(uPred_model(:,i) - u_l(:, i));
    Cost_model = Cost_model + Cost_model_vec(i);
end

Cost=0;
for i=1:maxsteps
    Cost_vec(i) = (uPred(:,i) - u_l(:, i))'*Q *(uPred(:,i) - u_l(:, i));
    Cost = Cost + Cost_vec(i);
end
meanZPCtime = mean(execTimeZPC)
stdZPCtime = std(execTimeZPC)
meanRMPCtime = mean(execTimeRMPC)
stdRMPCtime = std(execTimeRMPC)

%save the workspace
%save('workspaces\ZPC');
save('zonoDDSF\ddpc\workspaces\zonoDDSF');
%next run plotPolyZono for plotting