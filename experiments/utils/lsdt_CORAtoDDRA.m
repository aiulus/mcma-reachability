function [Y_full, U_full] = lsdt_CORAtoDDRA(testSuite)
    % convertTestSuiteToSequences - Converts testSuite into flat input and output matrices
    %
    % Inputs:
    %    testSuite - cell array of testCase objects
    %
    % Outputs:
    %    U_full - p x T matrix of control inputs
    %    Y_full - p x T matrix of outputs
    
    n_k = length(testSuite); % Number of unique input sequences
    assert(n_k > 0, 'testSuite must not be empty');
    
    % Initialize based on the first test case
    single_traj = testSuite{1};
    T_k = size(single_traj.y, 1); % time steps per trajectory
    n = size(single_traj.y, 2);
    m = size(single_traj.u, 2); % input dimension (assumed same as output for simplicity)
    T = n * n_k * T_k;
    
    U_full = zeros(m, T);
    Y_full = zeros(n, T);
    
    for i = 1:n_k
        tc = testSuite{i};
        %u_idx_start = (i - 1) * T_k + 1;
        %u_idx_end = i * T_k;

        %y_idx_start = (i - 1) * T_k + 1;
        %y_idx_end = i * T_k;

        lb = (i - 1) * T_k + 1;
        ub = i * T_k;
    
        % Ensure correct orientation: u and y are expected to be p x n_k
        u_i_1 = squeeze(tc.u(:, :, 1));
        u_i_2 = squeeze(tc.u(:, :, 2));
        y_i = squeeze(tc.y(:,:,1)); % assuming n_s = 1, remove singleton third dim
    
        %U_full(:, u_idx_start:u_idx_end) = u_i(:, :, i).';
        %Y_full(:, y_idx_start:y_idx_end) = y_i(:, :, i).';

        U_full(:, lb:ub) = u_i_1.';
        Y_full(:, lb:ub) = y_i.';
    end
end
