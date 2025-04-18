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
% Adapted by: Aybüke Ulusarslan, Yongkuan Zhang
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

%% initial set and input

r_u = 8; %reference input
%% TODO: add error handling for the inverse
r_y = u2y(r_u, sys); %reference output

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

%prepare M_Sigma which is a set of [A B]
M_Sigma = (Y_1T + -1* VmatZono + -1*WmatZono+VAmatZono)*pinv([Y_0T;U_data]);


%double check if the true A B is part of M_Sigma
intAB11 = intervalMatrix(M_Sigma);
intAB1 = intAB11.int;
intAB1.sup >= [sys_d.A,sys_d.B]
intAB1.inf <= [sys_d.A,sys_d.B]
% check the rank of the data
rank = rank([Y_0T;U_data])


%% Compute ZPC problem
%Horizon N for ZPC
N = 2;
%define output cost matrix
Qy = 1e3*eye(5); 
%control cost matrix
Qu = 0.001*eye(1);


execTimeZPC=[];
execTimeRMPC=[];
% ZPC number of time steps
maxsteps = 80;
% chosen time step for plotting 
chosedtimestep = 10;
for timesteps = 1:maxsteps
    if timesteps == 1
        % set the initial output to y0
        y_t(:,timesteps) = y0;
        y_t_model(:,timesteps) = y0;
        YPred(:,1) = y0;
    end
    
    
    % sdpvar variables
    u = sdpvar(1*ones(1,N),ones(1,N));
    y = sdpvar(5*ones(1,N+1),ones(1,N+1));
    alpha_u = sdpvar(1,N);
    sinf = sdpvar(5*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(5*ones(1,N+1),ones(1,N+1));
    R={};
    R{1} = zonotope([y_t(:,timesteps)]);
    %set the first constraint as y_t = current y
    Constraints = [y_t(:,timesteps) == y{1}];%,...
   
    
    for i = 1:N        
        genleni = size(R{i}.generators,2);
        %compute the reachable set for ZPC
        %card_cen = [R{i}.center;u{i}];        
        %card_zono = zonotope([card_cen,[R{i}.generators;zeros(1,genleni)]]);
        placeholder_u = r_u;  % safe numeric value
        dummy_cen = [R{i}.center; placeholder_u];
        dummy_gen = [R{i}.generators; zeros(1, genleni)];
        card_zono = zonotope(dummy_cen, dummy_gen);

        
        ABcard = [sys_d.A , sys_d.B]* card_zono;
        R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
       

        % AB * card_zono + W_sdp;
        %G_next = [ABcard.generators, W.generators, V.generators, AV.generators];
        % Use symbolic card_cen as center here
        %R{i+1} = zonotope(card_cen, G_next);

        % convert R to interval
        % extract center, determine left and right limit of the reahable set
        c = R{i+1}.c;
        delta = sum(abs(R{i+1}.G),2);
        
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit{i} = c - delta;
        rightLimit{i} = c + delta;
        
        %specify the constraints
        Constraints = [Constraints,...
            u{i} == U.center + alpha_u(i) * U.generators,...
            y{i+1} - sinf{i} == leftLimit{i},...
            y{i+1} + ssup{i} == rightLimit{i},...
            y{i+1}   - sinf{i} >= intc.inf,...
            y{i+1}   + ssup{i} <= intc.sup,...
            sinf{i} >= zeros(5,1),...
            ssup{i} >= zeros(5,1),...
            alpha_u(i) <= 1 , ...
            alpha_u(i) >= -1, ...
            ];
    end
    
    % chose the cost of ZPC
    Cost=0;
    for i=1:N
        Cost = Cost + (y{i+1}-r_y)'*Qy*(y{i+1}-r_y)+ (u{i}-r_u)'*Qu*(u{i}-r_u);
    end
    %solve ZPC
    options = sdpsettings('verbose',0,'solver','mosek');
    tic
    Problem = optimize(Constraints,Cost,options)
    execTimeZPC=[execTimeZPC,toc];
    Objective = double(Cost);
    uPred(timesteps) = double(u{1})
    YPred(:,timesteps+1) = double(y{2});
    %%

    %% save for plotting
    Rplotall{timesteps} = interval(zonotope(R{2}.c, R{2}.G));
    %%  ploting
    if chosedtimestep == timesteps
        for i =1:N+1
            RoverN{i} = zonotope(R{i}.c, R{i}.G);
            RoverN_int{i} = interval(RoverN{i});
            yoverN{i} =double(y{i});
            if i<N+1
                uoverN{i} =double(u{i});
            end
        end
    end
	
	
	
    %% ZPC given the model (RMPC-zono)
    % Control					
    alpha_u = sdpvar(1,N);
    sinf = sdpvar(5*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(5*ones(1,N+1),ones(1,N+1));
    R={};
    R{1} = zonotope([y_t_model(:,timesteps)]);
    u_model = sdpvar(1*ones(1,N),ones(1,N));
    y_model = sdpvar(5*ones(1,N+1),ones(1,N+1));
    Constraints = [y_t_model(:,timesteps) == y_model{1}];
    for i = 1:N
        %compute the reachable set for ZPC
        genleni = size(R{i}.generators,2);
        %card_cen = [R{i}.center;u{i}];        
        %card_zono = zonotope([card_cen,[R{i}.generators;zeros(1,genleni)]]);
        % Use numeric placeholder for center in zonotope creation
        placeholder_u = r_u;  % safe numeric value
        dummy_cen = [R{i}.center; placeholder_u];
        dummy_gen = [R{i}.generators; zeros(1, genleni)];
        card_zono = zonotope(dummy_cen, dummy_gen);

        ABcard = intervalMatrix(M_Sigma)* card_zono;
        R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
        %convert R to interval
        %extract center
        c = R{i+1}.Z(:,1);      
    
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit{i} = c - delta;
        rightLimit{i} = c + delta;
        
        
        Constraints = [Constraints,...
            u_model{i} == U.center + alpha_u(i) * U.generators,...
            y_model{i+1} - sinf{i} == leftLimit{i},...
            y_model{i+1} + ssup{i} == rightLimit{i},...
            y_model{i+1}   - sinf{i} >= intc.inf,...
            y_model{i+1}   + ssup{i} <= intc.sup,...
            sinf{i} >= zeros(5,1),...
            ssup{i} >= zeros(5,1),...
            alpha_u(i) <= 1 , ...
            alpha_u(i) >= -1, ...
            ];
    end
    
    
    
    Cost_model=0;
    for i=1:N
        Cost_model = Cost_model + (y_model{i+1}-r_y)'*Qy*(y_model{i+1}-r_y)+ (u_model{i}-r_u)'*Qu*(u_model{i}-r_u);
    end
    options = sdpsettings('verbose',0,'solver','mosek');
    tic
    Problem = optimize(Constraints,Cost_model,options)
    execTimeRMPC=[execTimeRMPC,toc];
    Objective = double(Cost_model);
    uPred_model(timesteps) = double(u_model{1});
    YPred_model(:,timesteps+1) = double(y_model{2});

    
    % apply the optimal control input to the plant 
    w_point = randPoint(W);
    v_point = randPoint(V);
    y_t(:,timesteps+1) = sys_d.A * y_t(:,timesteps) + sys_d.B * uPred(timesteps) + w_point +v_point - sys_d.A *v_point;
    y_t_model(:,timesteps+1) = sys_d.A * y_t_model(:,timesteps) + sys_d.B * uPred_model(timesteps) + w_point +v_point - sys_d.A *v_point;
    
    yt2ref(timesteps)= norm(y_t(:,timesteps)-r_y,2)
    yt2ref_model(timesteps)= norm(y_t_model(:,timesteps)-r_y,2)
    halt = 1;
end



Cost_model=0;
for i=1:timesteps
    Cost_model_vec(i) = (y_t_model(:,i+1)-r_y)'*Qy*(y_t_model(:,i+1)-r_y)+ (uPred_model(:,i)-r_u)'*Qu*(uPred_model(:,i)-r_u);
    Cost_model = Cost_model + Cost_model_vec(i);

end

Cost=0;
for i=1:timesteps
    Cost_vec(i) = (y_t(:,i+1)-r_y)'*Qy*(y_t(:,i+1)-r_y)+ (uPred(:,i)-r_u)'*Qu*(uPred(:,i)-r_u);
    Cost = Cost + Cost_vec(i);
end
meanZPCtime= mean(execTimeZPC)
stdZPCtime= std(execTimeZPC)
meanRMPCtime= mean(execTimeRMPC)
stdRMPCtime= std(execTimeRMPC)

%save the workspace
save('zonoDDSF\ddpc\workspaces\ZPC');
%next run plotPolyZono for plotting