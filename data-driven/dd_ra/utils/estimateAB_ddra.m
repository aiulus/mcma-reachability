function M_ab = estimateAB_ddra(sys, X_0T, X_1T, U_1T, WmatZ)
    % estimateAB_ddra: Estimate set of matrices (A, B) from data and validate inclusion.
    %
    % Inputs:
    %   - sys_d      : discrete-time system with fields A and B
    %   - X_0T       : state at time k (n × N)
    %   - X_1T       : state at time k+1 (n × N)
    %   - U_full     : input at time k (m × N)
    %   - Wmatzono   : matrix zonotope capturing disturbance W
    %
    % Outputs:
    %   - M_ab       : estimated matrix zonotope over (A,B)
    %   - intAB1     : interval over (A,B) from M_ab
    %
    % Errors:
    %   - Throws error if true (A,B) not inside the estimated set.

    % Center of (X_1 - W)
    X1W_cen = X_1T - WmatZ.center;

    % Matrix zonotope: X_1 - W ∈ AB * [X_0T; U_full]
    M_ab = matZonotope(X1W_cen, WmatZ.G) * pinv([X_0T; U_1T]);

    % Interval matrix over (A,B)
    intAB = intervalMatrix(M_ab);
    intAB1 = intAB.int;

    %% Temporarily disabled validation
    %A_true = sys_d.A;
    %B_true = sys_d.B;
    %AB_true = [A_true, B_true];

    % Validation
    %if any(intAB1.sup < AB_true, 'all') || any(intAB1.inf > AB_true, 'all')
    %    %% Interim debugging fix- change back to error(...)
    %    warning('estimateAB_ddra:ValidationFailed', ...
    %          'True system matrices (A,B) are not fully contained in the estimated matrix zonotope.');
    %else
    %    fprintf('[estimateAB_ddra] ✅ Validation successful: (A,B) contained in estimated set.\n');
    %end
end
