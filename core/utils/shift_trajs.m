function [X_0T, X_1T, U_0T, U_1T, U_full] = shift_trajs(testSuite)
% narx_CORAtoDDRA - Convert CORA testSuite to DDRA-compatible flat data format.
% Assumes y = x (full observability) and sys is a NARX system with .mFile handle

    %assert(isfield(sys, 'mFile') && isa(sys.mFile, 'function_handle'), ...
    %    'Expected sys to contain a function handle field "mFile".');

    K = numel(testSuite); % # Always equals one under the current createTestSuite()
    dim_y = size(testSuite{1}.y, 2);
    dim_u = size(testSuite{1}.u, 2);

    % s - #random samples per unique input sequence
    % T_k - trajectory length
    [T_k, ~, s] = size(testSuite{1}.u);
    T_k_minus = T_k - 1;
    totalTraj = K * s;

    x_all = zeros(dim_y, totalTraj);
    U_full = zeros(dim_u, totalTraj);

    X_0T = zeros(dim_y, totalTraj);
    X_1T = zeros(dim_y, totalTraj);
    U_0T = zeros(dim_u, totalTraj);
    U_1T = zeros(dim_u, totalTraj);

    for i = 1:K
        tc = testSuite{i};
        y_full = tc.y; % $\in (T_k \times n)$
        u_full = tc.u; % $\in (T_k \times m)$  

        for j = 1:totalTraj
            lb = (j - 1)*T_k + 1;
            ub = j*T_k;
            lb_minus = (j - 1)*T_k_minus + 1;
            ub_minus = j*T_k_minus;

            y_ij = squeeze(y_full(:, :, j)); % \in (T_k \times n)
            y_ij_plus = y_ij(2:end, :);
            y_ij_minus = y_ij(1:(end-1), :);
            x_all(:, lb:ub) = y_ij.';
            X_0T(:, lb_minus:ub_minus) = y_ij_minus.';
            X_1T(:, lb_minus:ub_minus) = y_ij_plus.';

            u_ij = squeeze(u_full(:, :, j));
            U_full(:, lb:ub) = u_ij.';
            if dim_u > 1
                u_ij_plus = u_ij(2:end, :);
                u_ij_minus = u_ij(1:(end-1), :);
            else
                u_ij_plus = zeros(T_k_minus, 1);
                u_ij_minus = zeros(T_k_minus, 1);
            end
            U_0T(:, lb_minus:ub_minus) = u_ij_minus.';
            U_1T(:, lb_minus:ub_minus) = u_ij_plus.';
        end
    end
    % Return transposed: T × n format
    %x_all = x_all.';
    %utraj_all = utraj_all.';
end

% ------------------------------ END OF CODE ------------------------------