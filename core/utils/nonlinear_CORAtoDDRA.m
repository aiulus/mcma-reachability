function [x_all, utraj_all] = nonlinear_CORAtoDDRA(testSuite)
% narx_CORAtoDDRA - Convert CORA testSuite to DDRA-compatible flat data format.
% Assumes y = x (full observability) and sys is a NARX system with .mFile handle

    %assert(isfield(sys, 'mFile') && isa(sys.mFile, 'function_handle'), ...
    %    'Expected sys to contain a function handle field "mFile".');

    n_m = numel(testSuite); % # Always equals one under the current createTestSuite()
    dim_y = size(testSuite{1}.y, 2);
    dim_u = size(testSuite{1}.u, 2);

    % s - #random samples per unique input sequence
    % T_k - trajectory length
    [n_k, ~, n_s] = size(testSuite{1}.y);

    x_all = zeros(dim_y, n_m * n_k * n_s);
    utraj_all = zeros(dim_u, n_m * n_k * n_s);

    for m = 1:n_m
        tc = testSuite{m};
        y_full = tc.y; % $\in (T_k \times n)$
        u_full = tc.u; % $\in (T_k \times m)$       

        for s = 1:n_s
            lb = (s-1)*n_k+1;
            ub = s*n_k;
            %u_ij = squeeze(u_full(:, :, j));
            u_ij = squeeze(u_full(:, :));
            if size(u_ij,1) ~= dim_u && size(u_ij,2) == dim_u
                u_ij = u_ij.';           
            end
            y_ij = squeeze(y_full(:, :, s)).'; % \in (T_k \times n)
            x_all(:, lb:ub) = y_ij;
            utraj_all(:, lb:ub) = u_ij;
        end
    end

    % Return transposed: T × n format
    %x_all = x_all.';
    %utraj_all = utraj_all.';
end
