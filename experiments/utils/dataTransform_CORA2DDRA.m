function [x_all, utraj_all] = dataTransform_CORA2DDRA(testSuite, n_s, n_m)
% narx_CORAtoDDRA - Convert CORA testSuite to DDRA-compatible flat data format.
% Assumes y = x (full observability) and sys is a NARX system with .mFile handle

    %assert(isfield(sys, 'mFile') && isa(sys.mFile, 'function_handle'), ...
    %    'Expected sys to contain a function handle field "mFile".');

    K = numel(testSuite); % # Always equals one under the current createTestSuite()
    dim_y = size(testSuite{1}.y, 2);

    % s - #random samples per unique input sequence
    [T_k, ~, s] = size(testSuite{1}.u);
    dim_u = size(testSuite{1}.u, 2);

    totalTraj = K * s;
    trajLength = T_k;

    x_all = zeros(dim_y, totalTraj * trajLength);
    utraj_all = zeros(dim_u, totalTraj * trajLength);

    trajIdx = 0;

    for i = 1:K
        u_full = tc.u;
        y_full = tc.y;
    
        if ndims(u_full) == 2
            u_full = repmat(u_full, [1 1 size(y_full,3)]);   % duplicate for every s
        end    

        for j = 1:totalTraj
            u_ij = squeeze(u_full(:, :, j));
            if size(u_ij,1) ~= dim_u && size(u_ij,2) == dim_u
                u_ij = u_ij.';                % make it (n_u × T_k)
            end
            trajIdx = trajIdx + 1;
            lb = (trajIdx-1)*T_k + 1;
            ub = trajIdx*T_k;
            x_all    (:, lb:ub) = squeeze(y_full(:,:,j)).';  
            utraj_all(:, lb:ub) = u_ij;      
        end
    end

    % Return transposed: T × n format
    %x_all = x_all.';
    %utraj_all = utraj_all.';
end


