function out = func_zonoDDSF(systype, dim, scale_noise, plot_toggle)
% Inputs:
%   system_type  - (char) system type name (e.g., 'chain_integrators')
%   dim          - (int) system dimensionality (only for scalable types)
%   scale_noise  - scaling for process noise
%   plot_toggle  - (logical) plot trajectory and reachable set projections

%------------- BEGIN CODE --------------

    %% INITIAL SETUP
    rng(4500);
    %clearvars;
    %close all;
    
    %% SYSTEM SETUP
    dt = 0.05; % sampling time
    sys = systemsZonoDDSF(systype, dt, dim);
    sys_c = sys.cont;
    sys_d = sys.discrete;
    n = sys.dims.n;
    
    %% SIMULATION SETUP
    initpoints =100; % number of trajectories
    steps =5; % number of steps for each trajectory
    totalsamples = initpoints*steps; %Total number of samples
    
    %% RETRIEVE BOUNDARY CONDITIONS
    sys = setupBoundaryConditions(sys);
    %intc = sys.bcs.intc;
    y0 = sys.bcs.y0;
    X0 = sys.bcs.X0;
    U = sys.bcs.U;
    
    %% SET UP NOISE TERMS
    w_spread = 0.01 * scale_noise.pcs; % controlls process noise levels
    v_spread = 0.002 * scale_noise.msmt; % controlls measurement noise levels
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
    if plot_toggle
        figure;
        subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
        subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
        close;
    end
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
    maxsteps = 10;
    % time step for plotting 
    timestep_plot = 10;
    
    [uPred, uPred_model, u_l_hist, y_t, y_t_model, execTimeZPC, execTimeRMPC] =  runZonoDDSF( ...
        sys, y0, U, N, maxsteps, timestep_plot, Q, W, V, AV, M_Sigma);          
    
    Cost_model_vec = sum((uPred_model - u_l_hist).^2 .* Q, 1);
    Cost_model = sum(Cost_model_vec);

    Cost_vec = sum((uPred - u_l_hist).^2 .* Q, 1);
    Cost = sum(Cost_vec);

    if size(y_t, 1) >=2
        lb = sys.bcs.intc.inf(2);
        ub = sys.bcs.intc.sup(2);
        c = conservatism(y_t(2, :), y_t_model(2, :), lb, ub);
    else
        c = NaN;
    end

    meanZPCtime = mean(execTimeZPC)
    stdZPCtime = std(execTimeZPC)
    meanRMPCtime = mean(execTimeRMPC)
    stdRMPCtime = std(execTimeRMPC)
    
    %save the workspace
    %save('workspaces\ZPC');
    outputFolder = fullfile('zonoDDSF', 'workspaces', 'data-driven', 'pc');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    save('zonoDDSF\workspaces\data-driven\pc\zonoDDSF');

    out = struct( ...
        'Cost', Cost, ...
        'Cost_model', Cost_model, ...
        'execTimeZPC', execTimeZPC, ...
        'execTimeRMPC', execTimeRMPC, ...
        'conservatism', c, ...
        'y_t', y_t, ...
        'y_t_model', y_t_model, ...
        'uPred', uPred, ...
        'uPred_model', uPred_model ...
    );

    %next run plotPolyZono/plotZonoDDSF for plotting
    %plotZonoDDSF;
    %plotPolyZono;
end
