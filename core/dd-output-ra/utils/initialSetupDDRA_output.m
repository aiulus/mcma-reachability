function [X0, U, W, M_w, V, M_v, AV, M_Av] = initialSetupDDRA_output(sys, initpoints, steps,  X0_center, X0_spread, ...
    U_center, U_spread, W_center, W_spread, V_center, V_spread, x0GenOrder, uGenOrder)
    
    totalsamples = initpoints * steps;
    X0 = zonotope(X0_center, X0_spread*ones(sys.nrOfStates, x0GenOrder));
    U = zonotope(U_center, U_spread*ones(sys.nrOfInputs, uGenOrder));
    % Achtung-- Process noise generator order hard-coded to 1
    W = zonotope(W_center, W_spread*ones(sys.nrOfStates, 1));
    for i=1:size(W.gen, 2)
        vec = W.G(:, i);
        GW{i} = [vec, zeros(sys.nrOfStates, totalsamples)];
        for j=1:totalsamples
            GW{j+1} = [GW{i+j-1}(:, 2:end) GW{i+j-1}(:, 1)];
        end
    end
    M_w = matZonotope(zeros(sys.nrOfStates, totalsamples), GW);
    % Achtung-- Measurement noise generator order hard-coded to 1
    V = zonotope(V_center, V_spread*ones(sys.nrOfStates, 1));
    CV = zeros(sys.nrOfStates, totalsamples);
    
    for i=1:size(V.G, 2)
        vec = vec.G(:, i);
        GV{i} = [vec, zeros(sys.nrOfStates, totalsamples - 1)];
        for j=1:totalsamples-1
            GV{j+i} = [GV{i+j-1}(:, 2:end) GV{i+j-1}(:, 1)];
        end
    end
    M_v = matZonotope(CV, GV);
    
    AV = sys_d.A*V;
    M_Av = sys_d.A*M_v;
end

