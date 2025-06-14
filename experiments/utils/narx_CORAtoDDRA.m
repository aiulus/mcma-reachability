function [x_all, utraj_all] = aux_CORAtoDDRA(testSuite, sys)
% AUX_CORATODDRA Convert CORA testSuite to DDRA flat data for nonlinearARX
% Assumes y = x (full observability) and sys is nonlinearARX

    N = numel(testSuite);
    dim_y = sys.nrOfOutputs;
    dim_u = sys.nrOfInputs;
    n_p = sys.n_p;
    [T, ~, s] = size(testSuite{1}.u);

    totalTraj = N * s;
    x_all     = zeros(totalTraj * dim_y, T+1);
    utraj_all = zeros(totalTraj * dim_u, T);

    for i = 1:N
        tc = testSuite{i};
        Y0 = tc.initialState;   % dim_y×1×s
        U = tc.u;               % T×dim_u×s

        for j = 1:s
            traj_idx = (i-1)*s + j;
            rows_x = (traj_idx-1)*dim_y + (1:dim_y);
            rows_u = (traj_idx-1)*dim_u + (1:dim_u);

            y_hist = zeros(dim_y, n_p);
            x_traj = zeros(dim_y, T+1);
            u_traj = squeeze(U(:,:,j))'; % dim_u×T
            x_traj(:,1) = Y0(:,:,j);

            % Fill initial history with duplicates of x0 if no prior data
            for p = 1:n_p
                y_hist(:,p) = Y0(:,:,j);  % or use zeros if preferred
            end

            for k = 1:T
                u_hist = zeros(dim_u, n_p+1);
                for p = 0:n_p
                    if (k - p) >= 1
                        u_hist(:,p+1) = u_traj(:,k - p);
                    else
                        u_hist(:,p+1) = u_traj(:,1); % replicate first input
                    end
                end

                y_input = y_hist;
                u_input = u_hist;

                % Evaluate NARX: y(k) = f(y(k-1...k-n_p), u(k)...u(k-n_p))
                y_new = sys.mFile(y_input, u_input);

                x_traj(:,k+1) = y_new;

                % Update history
                y_hist = [y_hist(:,2:end), y_new];
            end

            x_all(rows_x, :)     = x_traj;
            utraj_all(rows_u, :) = u_traj;
        end
    end
end
