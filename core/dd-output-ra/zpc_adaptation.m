% Adapted from:
% https://github.com/aalanwar/Data-Driven-Predictive-Control/ZPC.m
rand('seed', 1);

clear all;
%close all;

%% ------------- 1.1) Define the system -------------
dim_x = 5;

% System in cont time
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
C = [1,0,0,0,0];
D = 0;
sys_c = ss(A,B_ss,C,D);

% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);

% extract dimensions
n = size(sys_d.A, 1); m = size(sys_d.B, 2);

%% ------------- 1.2) Data properties -------------
initpoints = 100; % #(trajectories), maps to n_m
steps = 5; % traj. length, maps to n_k
totalsamples = initpoints*steps;

%% ------------- 1.3) Problem setting -------------
uref = 8; % reference input
ref = inv(eye(5)-sys_d.A)*sys_d.B*uref; % reference output

% output constraintS
y_lb = [-10;2;-10;-10;-10]; 
y_ub = [10;10;10;10;10]; 
intc = interval(y_lb,y_ub);

%initial point
y0 = [-2;4;3;-2.5;5.5];

%% Bounded uncertainty models
[X0, U, W, M_w, V, M_v, AV, M_Av] = initialSetupDDRA_output(sys_d, initpoints, steps,  y0, 25, ...
    uref-1, 20-1, zeros(n, 1), 0.01, zeros(n, 1), 0.002, 1, 1);

%% Control input sequence
for i=1:totalsamples
    u_seq(i) = randPoint(U);
end

[ytraj, ytraj_v, utraj] = getDataDDRA_output(sys_d, initpoints, steps, X0, W, V, u_seq);

U_data = utraj(:,1:totalsamples); 
Y_0T = ytraj_v(:,1:totalsamples);
Y_1T = ytraj_v(:,1:totalsamples);

M_sigma = estimateAB_output(sys_d, Y_0T, Y_1T, U_0T, M_v, M_w, M_Av);

[R_tr, R_pred] = propDDRA_output(sys_d, u_seq, X0, W, V, AV, M_sigma, timesteps);