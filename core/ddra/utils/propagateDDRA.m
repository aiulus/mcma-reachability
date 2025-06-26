function [X_model, X_data] = propagateDDRA(U_full, X_0T, X_1T, X0, U, W, sys, M_ab, totalsteps)
    % propagateDDRA: Propagates reachable sets using model-based and data-driven methods.
    %
    % Inputs:
    %   - X0         : initial set (zonotope)
    %   - U          : input set (zonotope)
    %   - W          : disturbance set (zonotope)
    %   - sys_d      : discrete-time system with fields A, B
    %   - M_ab       : matrix zonotope estimate of [A B]
    %   - totalsteps : number of propagation steps
    %
    % Outputs:
    %   - X_model    : cell array of reachable sets using the known model
    %   - X_data     : cell array of reachable sets using data-driven estimate
    if isa(sys, 'linearSysDT') || isa(sys, 'linearSys')
        % Initialize reachable set containers
        X_model = cell(totalsteps + 1, 1);
        X_data = cell(totalsteps + 1, 1);
    
        % Set initial condition
        X_model{1} = X0;
        X_data{1} = X0;
    
        for i = 1:totalsteps
            fprintf('Computing step %d / %d of reachset propagation...\n', i, totalsteps);
            % Reduce for computational efficiency (model-based)
            X_model{i} = reduce(X_model{i}, 'girard', 400);             
            X_model{i+1} = sys.A * X_model{i} + sys.B * U + W;
            % Reduce for computational efficiency (data-driven)
            X_data{i} = reduce(X_data{i}, 'girard', 400);
            X_data{i+1} = M_ab * cartProd(X_data{i}, U) + W;
        end
    elseif isa(sys, 'nonlinearARX') || isa(sys, 'linearARX')
        cfg = getConfig();
        options_reach = cfg.options_reach;
        %% Used in linearize_DT (line 34), not sure what it does
        % params = struct('R0', X0, 'U', U, 'tFinal', totalsteps);
        options = struct( ...
        'dim_x', sys.nrOfStates, ...
        'U_full', U_full, ...
        'X_0T', X_0T, ...
        'X_1T', X_1T, ...
        'W', W, ...
        'R0', X0, ...
        'U', U, ...
        'tFinal', totalsteps, ...
        'zonotopeOrder', options_reach.zonotopeOrder, ...
        'tensorOrder', options_reach.tensorOrder, ...
        'errorOrder', options_reach.errorOrder ...
        );
        
        %% Achtung! sys.B_bar{2} hardcoded for n_p = 1
        %X_model{i+1} = sys.A_bar * X_model{i} + sys.B_bar{2} * U + W;     
        %X_model = reach_DT(sys, params, options);  
        %X_data = cell(totalsteps + 1, 1);        
        %X_data{1} = X0; % Set initial condition

        % Initialize reachable set containers
        X_model = cell(totalsteps + 1, 1);
        X_data = cell(totalsteps + 1, 1);
    
        % Set initial condition
        X_model{1} = X0;
        X_data{1} = X0;
    
        for i = 1:totalsteps
            fprintf('Computing step %d / %d of reachset propagation...\n', i, totalsteps);
            % Reduce for computational efficiency (data-driven)
            % X_model{i + 1} = linReach_DT(sys, X_model{i}, options); 
            X_model{i+1} = cell2mat(sys.A_bar) * X_model{i} + sys.B_bar{2} * U + W;  
            X_data{i} = reduce(X_data{i}, 'girard', 400);
            X_data{i+1} = M_ab * cartProd(X_data{i}, U) + W;
        end
    end    
end
