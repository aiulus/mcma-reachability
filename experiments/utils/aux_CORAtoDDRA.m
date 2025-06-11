function [x_all, utraj_all] = aux_CORAtoDDRA(testSuite, sys)
% AUX_CORATODDRA  Convert CORA testSuite to flat data for DDRA pipeline
%   Supports multiple replicates per testCase by flattening along third dim.
% Inputs:
%   testSuite - 1×N cell of testCase objects, each with fields:
%                  .initialState  (n×1×s)
%                  .u             (T×m×s)
%   sys       - system struct with sys.discrete.A (n×n) and B (n×m)
% Outputs:
%   x_all     - (N*s*n)×(T+1) matrix: stacked state trajectories
%   utraj_all - (N*s*m)×T       matrix: stacked input trajectories

    N = numel(testSuite);
    % State dimension
    n = size(testSuite{1}.initialState,1);
    % Input trajectory size: T×m×s
    [T, m, s] = size(testSuite{1}.u);

    % Total trajectories = N * s
    totalTraj = N * s;

    % Preallocate
    x_all     = zeros(totalTraj * n, T+1);
    utraj_all = zeros(totalTraj * m, T);

    idx = 0;
    for i = 1:N
        tc = testSuite{i};
        % initialState: n×1×s
        X0 = tc.initialState;
        % u: T×m×s
        Uall = tc.u;
        for j = 1:s
            idx = idx + 1;
            x0 = X0(:,:,j);         % n×1
            Umat = squeeze(Uall(:,:,j));  % T×m
            % Reconstruct state trajectory
            Xmat = zeros(n, T+1);
            Xmat(:,1) = x0;
            for k = 1:T
                Xmat(:,k+1) = sys.A * Xmat(:,k) + sys.B * Umat(k,:)';
            end
            % Store
            rows_x = (idx-1)*n + (1:n);
            rows_u = (idx-1)*m + (1:m);
            x_all(rows_x, :)     = Xmat;
            utraj_all(rows_u, :) = Umat';
        end
    end
end
