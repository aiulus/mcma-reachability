function [X_0T, X_1T, U_0T, U_1T] = shift_trajs(testSuite)
% narx_CORAtoDDRA - Convert CORA testSuite to DDRA-compatible flat data format.
% Assumes y = x (full observability) and sys is a NARX system with .mFile handle

    %assert(isfield(sys, 'mFile') && isa(sys.mFile, 'function_handle'), ...
    %    'Expected sys to contain a function handle field "mFile".');

    K = numel(testSuite); % # distinct input sequences in the ensemble
    dim_y = size(testSuite{1}.y, 2);
    dim_u = size(testSuite{1}.u, 2);

    % s - #random samples per unique input sequence
    % T_k - trajectory length
    [T_k, ~, s] = size(testSuite{1}.u);
    totalTraj = K * s;

    x_all = zeros(dim_y, totalTraj);
    utraj_all = zeros(dim_u, totalTraj);

    X_0T = zeros(dim_y, totalTraj);
    X_1T = zeros(dim_y, totalTraj);
    U_0T = zeros(dim_u, totalTraj);
    U_1T = zeros(dim_u, totalTraj);

    for i = 1:K
        tc = testSuite{i};
        y_i = tc.y; % $\in (T_k \times n)$
        u_i = tc.u.'; % $\in (T_k \times m)$  
        u_i_plus = u_i(:, 2:end);
        u_i_minus = u_i(:, 1:(end-1));

        for j = 1:s
            lb = (s+i-1)*T_k+1;
            ub = (s+i)*T_k;
            lb_minus_one = (s+i-1)*(T_k - 1)+1;
            ub_minus_one = (s+i)*(T_k - 1);
            utraj_all(:, lb:ub) = u_i;
            y_ij = squeeze(y_i(:, :, j)).'; % \in (T_k \times n)
            y_ij_plus = y_ij(:, 2:end);
            y_ij_minus = y_ij(:, 1:(end-1));
            x_all(:, lb:ub) = y_ij;
            X_0T(:, lb_minus_one:ub_minus_one) = y_ij_minus;
            X_1T(:, lb_minus_one:ub_minus_one) = y_ij_plus;
            U_0T(:, lb_minus_one:ub_minus_one) = u_i_minus;
            U_1T(:, lb_minus_one:ub_minus_one) = u_i_plus;
        end
    end

    % Return transposed: T × n format
    x_all = x_all.';
    utraj_all = utraj_all.';
end
