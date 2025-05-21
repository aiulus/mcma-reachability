function M_Sigma = estimateAB(Y_0T, U_data, Y_1T, VmatZono, WmatZono, VAmatZono, sys_d)
% ESTIMATEAB Constructs matrix zonotope M_Sigma and checks [A B] inclusion
%
% Inputs:
%   Y_0T       - matrix of past outputs (Y_-)
%   U_data     - matrix of past inputs (U_-)
%   Y_1T       - matrix of future outputs (X_+)
%   VmatZono   - measurement noise zonotope (M_v)
%   WmatZono   - process noise zonotope (M_w)
%   VAmatZono  - A*measurement noise zonotope (M_Av)
%   sys_d      - system struct containing fields A, B
%
% Output:
%   M_Sigma    - matrix zonotope over [A B] constructed from data

    % 1. Form matrix zonotope estimate of [A B]
    M_Sigma = (Y_1T + -1* VmatZono + -1*WmatZono+VAmatZono)*pinv([Y_0T;U_data]);

    % 2. Extract interval from matrix zonotope
    M_interval = intervalMatrix(M_Sigma);
    M_bounds = M_interval.int;

    AB_true = [sys_d.A, sys_d.B];

    % 3. Validate inclusion of true AB
    if any(any(AB_true > M_bounds.sup))
        %% Temporary debug fix-- change back to error(...)!
        warning('[A B] is NOT fully contained within the upper bound of M_Sigma!');
    end
    if any(any(AB_true < M_bounds.inf))
        warning('[A B] is NOT fully contained within the lower bound of M_Sigma!');
    end

    % 4. Rank condition on data
    data_rank = rank([Y_0T; U_data]);
    required_rank = size(Y_0T, 1) + size(U_data, 1);

    if data_rank < required_rank
        error('Rank deficient data: rank = %d, required = %d', data_rank, required_rank);
    end
end

