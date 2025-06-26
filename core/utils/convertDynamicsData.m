function [x_all, utraj_all] = convertDynamicsData(testSuite, sys)
% CONVERTDYNAMICSDATA  Convert createTestSuite output into DDRA flat data
%   [x_all, utraj_all] = convertDynamicsData(testSuite, sys)
%
% Inputs:
%   testSuite - 1×N cell array of testCase objects (from createTestSuite)
%   sys       - system object returned by loadDynamics (must have .discrete.A,B)
%
% Outputs:
%   x_all     - (N*n)×(T+1) matrix of states for DDRA
%   utraj_all - (N*m)×T       matrix of inputs for DDRA
%
% Example:
%   [sys, R0, U, p_true] = loadDynamics('Square','rand');
%   ts = createTestSuite(sys, params, n_k, n_m, n_s, options);
%   [x_all, utraj_all] = convertDynamicsData(ts, sys);

    N = numel(testSuite);
    % Determine state/input dims and horizon
    n = size(testSuite{1}.initialState,1);
    if isempty(testSuite{1}.u)
        error('convertDynamicsData:NoInput','Each testSuite entry must have .u');
    end
    % testSuite{1}.u is T×m×replicates (rep=1)
    [T, m, rep] = size(testSuite{1}.u);
    if rep~=1
        error('convertDynamicsData:MultiRep','Only single replicate supported');
    end

    % Preallocate
    x_all     = zeros(N*n, T+1);
    utraj_all = zeros(N*m, T);

    for i = 1:N
        tc = testSuite{i};
        % initial state
        x0 = tc.initialState;
        % input trajectory: T×m
        Umat = squeeze(tc.u);
        % reconstruct state trajectory via discrete A,B
        Xmat = zeros(n, T+1);
        Xmat(:,1) = x0;
        for k = 1:T
            Xmat(:,k+1) = sys.discrete.A * Xmat(:,k) + sys.discrete.B * Umat(k,:)';
        end
        % store into flat arrays
        rows_x = (i-1)*n + (1:n);
        rows_u = (i-1)*m + (1:m);
        x_all(rows_x, :)     = Xmat;
        utraj_all(rows_u, :) = Umat';
    end
end
