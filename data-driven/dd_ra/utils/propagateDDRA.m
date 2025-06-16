function [X_model, X_data] = propagateDDRA(X0, U, W, sys, M_ab, totalsteps)
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
        if isa(sys, 'nonlinearARX')
            X_model{i+1} = sys.mFile(X_model{i}, U) + W;
        else
            X_model{i+1} = sys.A * X_model{i} + sys.B * U + W;
        end

        % Reduce for computational efficiency (data-driven)
        X_data{i} = reduce(X_data{i}, 'girard', 400);
        X_data{i+1} = M_ab * cartProd(X_data{i}, U) + W;
    end

end
