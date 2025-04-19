function [u_opt, y_next, u, y, execTime, R, R_plot] = ...
    solveQPzonoDDSF(R, current_y, sys_d, N, r_u, r_y, U, Qy, Qu, W, V, AV, intc, k, maxsteps, n)
    % Inputs match the ZPC planning section in the original function
    % Outputs:
    %   u_opt   - optimal control at step k (scalar)
    %   y_next  - predicted output y(:,k+1)
    %   execTime - time to solve optimization
    %   R_plot - interval zonotope for plotting (R{2})

    % ---- sdpvar definitions (copied from original ZPC section) ----
    u       = sdpvar(1*ones(1,N), ones(1,N));
    y       = sdpvar(n*ones(1,N+1), ones(1,N+1));
    alpha_u = sdpvar(1,N);
    s_inf    = sdpvar(n*ones(1,N+1), ones(1,N+1));
    s_sup    = sdpvar(n*ones(1,N+1), ones(1,N+1));

    Constraints = current_y == y{1};      % Original: Constraints = y_t(:,k) == y{1};

    leftLimit  = cell(1, N);
    rightLimit = cell(1, N);

    for i = 1:N
        genleni = size(R{i}.generators, 2);   % Original: genleni = size(R{i}.generators,2);
        dummy_cen = [R{i}.center; r_u];       % r_u is placeholder
        dummy_gen = [R{i}.generators; zeros(1, genleni)];
        card_zono = zonotope(dummy_cen, dummy_gen);    % Original: card_zono = zonotope(...);

        ABcard = [sys_d.A, sys_d.B] * card_zono;       % Original: ABcard = [sys_d.A , sys_d.B]*card_zono;
        R{i+1} = zonotope([ABcard.center, ...
            [ABcard.generators, W.generators, V.generators, AV.generators]]);

        c = R{i+1}.c;                                   % Original: c = R{i+1}.c;
        delta = sum(abs(R{i+1}.Z), 2) - abs(c);         % delta computation
        leftLimit{i}  = c - delta;
        rightLimit{i} = c + delta;

        % ---- Constraints as in original ----
        Constraints = [Constraints, ...
            u{i} == U.center + alpha_u(i) * U.generators, ...
            y{i+1} - s_inf{i} == leftLimit{i}, ...
            y{i+1} + s_sup{i} == rightLimit{i}, ...
            y{i+1} - s_inf{i} >= intc.inf, ...
            y{i+1} + s_sup{i} <= intc.sup, ...
            s_inf{i} >= zeros(n,1), ...
            s_sup{i} >= zeros(n,1), ...
            alpha_u(i) <= 1, ...
            alpha_u(i) >= -1];
    end

    %% Q_u, r_u, r_y
    % ---- Cost function ----
    Cost = 0;
    for i = 1:N
        Cost = (u{i} - u_l)' * Q * (u{i} - u_l);
    end

    % ---- Solve the optimization ----
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    tic;

    Problem = optimize(Constraints, Cost, options);
    fprintf("Problem (step %d / %d): \n", k, maxsteps);
    disp(Problem);

    Objective = double(Cost);
    fprintf("Objective (step %d / %d): \n", k, maxsteps);
    disp(Objective);
    
    execTime = toc;

    if Problem.problem ~= 0
        warning("ZPC optimization problem at step failed: %s", yalmiperror(Problem.problem));
    end

    % ---- Outputs (original lines) ----
    u_opt  = double(u{1});          % Original: uPred(k) = double(u{1});
    y_next = double(y{2});          % Original: YPred(:,k+1) = double(y{2});
    R_plot = interval(zonotope(R{2}.c, R{2}.G));  % Original: Rplotall{k} = ...
end


