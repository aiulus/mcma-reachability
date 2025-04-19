function [u_opt_model, y_next_model, execTime] = solveRMPC(current_y_model, M_Sigma, N, r_u, r_y, U, Qy, Qu, W, V, AV, intc, n)
    % Inputs correspond to the RMPC block in original function
    % Outputs:
    %   u_opt_model - optimal control at step k
    %   y_next_model - predicted output for model system
    %   execTime - solve time for optimization

    % ---- sdpvar variables (copied from RMPC block) ----
    alpha_u = sdpvar(1, N);
    sinf    = sdpvar(n*ones(1,N+1), ones(1,N+1));
    ssup    = sdpvar(n*ones(1,N+1), ones(1,N+1));
    u_model = sdpvar(1*ones(1,N), ones(1,N));
    y_model = sdpvar(n*ones(1,N+1), ones(1,N+1));
    % ---- Initial reachable set ----
    R = cell(1, N+1);
    R{1} = zonotope([current_y_model]);           % Original: R{1} = zonotope([y_t_model(:,k)]);
    Constraints = current_y_model == y_model{1};   % Original: Constraints = y_t_model(:,k) == y_model{1};

    leftLimit  = cell(1, N);
    rightLimit = cell(1, N);

    for i = 1:N
        genleni = size(R{i}.generators, 2);
        dummy_cen = [R{i}.center; r_u];
        dummy_gen = [R{i}.generators; zeros(1, genleni)];
        card_zono = zonotope(dummy_cen, dummy_gen);

        ABcard = intervalMatrix(M_Sigma) * card_zono;  % <-- Difference from ZPC
        R{i+1} = zonotope([ABcard.center, ...
            [ABcard.generators, W.generators, V.generators, AV.generators]]);

        c = R{i+1}.Z(:, 1);       % Original: c = R{i+1}.Z(:,1);
        delta = sum(abs(R{i+1}.Z), 2) - abs(c);
        leftLimit{i}  = c - delta;
        rightLimit{i} = c + delta;

        % ---- Constraints ----
        Constraints = [Constraints, ...
            u_model{i} == U.center + alpha_u(i) * U.generators, ...
            y_model{i+1} - sinf{i} == leftLimit{i}, ...
            y_model{i+1} + ssup{i} == rightLimit{i}, ...
            y_model{i+1} - sinf{i} >= intc.inf, ...
            y_model{i+1} + ssup{i} <= intc.sup, ...
            sinf{i} >= zeros(n,1), ...
            ssup{i} >= zeros(n,1), ...
            alpha_u(i) <= 1, ...
            alpha_u(i) >= -1];
    end

    % ---- Cost function ----
    Cost_model = 0;
    for i = 1:N
        Cost_model = Cost_model + ...
            (y_model{i+1} - r_y)' * Qy * (y_model{i+1} - r_y) + ...
            (u_model{i} - r_u)' * Qu * (u_model{i} - r_u);
    end

    % ---- Solve optimization ----
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    tic;
    Problem = optimize(Constraints, Cost_model, options);
    execTime = toc;

    if Problem.problem ~= 0
        warning("RMPC optimization problem failed: %s", yalmiperror(Problem.problem));
    end

    u_opt_model  = double(u_model{1});       % Original: uPred_model(k) = ...
    y_next_model = double(y_model{2});       % Original: YPred_model(:,k+1) = ...
end
