function [x_all, utraj_all] = narx_CORAtoDDRA(testSuite)
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

    for i = 1:K
        tc = testSuite{i};
        y_i = tc.y; % $\in (T_k \times n)$
        u_i = tc.u.'; % $\in (T_k \times m)$       

        for j = 1:s
            lb = (s+i-1)*T_k+1;
            ub = (s+i)*T_k;
            utraj_all(:, lb:ub) = u_i;
            y_ij = squeeze(y_i(:, :, j)).'; % \in (T_k \times n)
            x_all(:, lb:ub) = y_ij;
        end
    end

    % Return transposed: T × n format
    x_all = x_all.';
    utraj_all = utraj_all.';
end
