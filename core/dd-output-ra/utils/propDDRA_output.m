function [R_tr, R_pred] = propDDRA_output(sys, utraj, X0, W, V, AV, M_sigma, timesteps)
    dim_u = sys.nrOfInputs;
    dim_y = sys.nrOfOutputs;

    % y0 is hard-coded and used as the center vector to build X0 in
    % the original code (Data-Driven-Predictive-Control\ZPC.m)
    y0 = X0.c;

    % Initialize containers
    R_tr = cell(timesteps + 1, 1);
    R_pred = cell(timesteps + 1, 1);
    R_pred{1} = zonotope(y0);
    R_tr{1} = zonotope(y0);

    for t=1:timesteps
        card_cen = [R_pred{t}.c; utraj{t}];
        genVec_count = size(R_pred{t}.G, 2);
        card_zono = zonotope(card_cen, [R_pred{t}.G; zeros(dim_u, genVec_count)]); % is 1 hard-coded for dim_u?
        ABcard = intervalMatrix(M_sigma) * card_zono;
        R_pred{t+1} = zonotope(ABcard.c, [ABcard.G, W.G, V.G, AV.G]);        
    end

    for t=1:timesteps
        card_cen = [R_tr{t}.c; utraj{t}];
        genVec_count = size(R_tr{t}.G, 2);
        card_zono = zonotope(card_cen, [R_tr{t}.G; zeros(dim_u, genVec_count)]); % is 1 hard-coded for dim_u?
        ABcard = [sys.A, sys.B] * card_zono;
        R_tr{t+1} = zonotope(ABcard.c, [ABcard.G, W.G, V.G, AV.G]);        
    end
end

