function [X0, U, W, M_w, V, M_v, AV, M_Av] = initialSetupDDRA_output( ...
    sys, initpoints, steps, n_s, ...
    X0_center, X0_spread, ...
    U_center, U_spread, ...
    W_center, W_spread, ...
    V_center, V_spread, ...
    x0GenOrder, uGenOrder)

    % System dimensions, assuming p = n (CORA)
    %n = sys.nrOfStates;
    %m = sys.nrOfInputs;

    % System dimensions, assuming p = n (ss)
    n = size(sys.A, 1);
    m = size(sys.B, 2);
    totalsamples = initpoints * (steps - 1);

    % Initial output set
    X0 = zonotope(X0_center, X0_spread * ones(n, x0GenOrder));

    % Input set
    U = zonotope(U_center, U_spread * ones(m, uGenOrder));

    % Process noise W
    W = zonotope(W_center, W_spread * ones(n, 1));

    for i = 1:size(W.G, 2)
        vec = W.G(:, i);
        GW{i} = [vec, zeros(n, totalsamples - 1)];
        for j = 1:totalsamples - 1
            GW{j + i} = [GW{i + j - 1}(:, 2:end), GW{i + j - 1}(:, 1)];
        end
    end
    M_w = matZonotope(zeros(n, totalsamples), GW);

    % Measurement noise V
    V = zonotope(V_center, V_spread * ones(n, 1));
    CV = zeros(n, totalsamples);

    for i = 1:size(V.G, 2)
        vec = V.G(:, i);
        GV{i} = [vec, zeros(n, totalsamples - 1)];
        for j = 1:totalsamples - 1
            GV{j + i} = [GV{i + j - 1}(:, 2:end), GV{i + j - 1}(:, 1)];
        end
    end
    M_v = matZonotope(CV, GV);

    % Compute A*V and A*M_v
    AV = sys.A * V;
    M_Av = sys.A * M_v;
end
