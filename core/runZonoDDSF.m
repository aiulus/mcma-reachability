function [uPred, uPred_model, y_t, y_t_model, execTimeZPC, execTimeRMPC] = ...
    runZonoDDSF(sys_d, y0, intc, U, N, maxsteps, timestep_plot, Q, W, V, AV, M_Sigma)
    
    %% PREALLOCATION STEPS
    p = size(y0,1);
    y_t = zeros(p, maxsteps + 1);
    y_t_model = zeros(p, maxsteps + 1);

    uPred        = zeros(1, maxsteps);
    uPred_model  = zeros(1, maxsteps);

    YPred        = zeros(p, maxsteps + 1);
    YPred_model  = zeros(p, maxsteps + 1);

    execTimeZPC  = zeros(1, maxsteps);
    execTimeRMPC = zeros(1, maxsteps);

    Rplotall     = cell(1, maxsteps);

    yt2ref       = zeros(1, maxsteps);
    yt2ref_model = zeros(1, maxsteps);
    
    % For optional plotting at timestep_plot
    R_N = cell(1, N+1);
    R_N_int = cell(1, N+1);
    y_N = cell(1, N+1);
    u_N = cell(1, N);

    for k = 1:maxsteps
        if k == 1
            [y_t, y_t_model, YPred] = runZPCinitializeY0(y0, y_t, y_t_model, YPred, k);
        end

        % ---- Initialize reachable sets ----
        R = cell(1, N+1);
        R{1} = zonotope(y_t(:,k));         % Original: R{1} = zonotope([y_t(:,k)]);

        %% Compute ZPC problem
        [uPred(k), YPred(:,k+1), u, y, execTimeZPC(k), R, Rplotall{k}] = ...
            solveZDDSF(R, y_t(:,k), sys_d, N, Q, W, V, AV, intc, k, maxsteps, p);

        %%  Plot
        if timestep_plot == k
            for i =1:N+1
                R_N{i} = zonotope(R{i}.c, R{i}.G);
                R_N_int{i} = interval(R_N{i});
                y_N{i} =double(y{i});
                if i<N+1
                    u_N{i} =double(u{i});
                end
            end
        end  

        %% ZPC given the model (RMPC-zono)
        [uPred_model(k), YPred_model(:,k+1), execTimeRMPC(k)] = ...
             solveRDDSF(y_t_model(:,k), M_Sigma, N,  U, Q, W, V, AV, intc, p);	
        
        %% Apply the optimal control input
        [y_t(:,k+1), y_t_model(:,k+1)] = ...
            sysIter(y_t(:,k), y_t_model(:,k), uPred(k), uPred_model(k), sys_d, W, V);
        
        yt2ref(k) = norm(y_t(:,k)-r_y,2);
        yt2ref_model(k) = norm(y_t_model(:,k)-r_y,2);

        %% TODO: Find out why this is here
        % halt = 1;
    end
end

