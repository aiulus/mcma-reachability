function [U_data, Y_0T, Y_1T] = prepareDataMatrices(x, x_v, utraj, dim_x, steps, total_samples)
% prepareDataMatrices - Generates Hankel-style training data
%
% Inputs:
%   x        - Clean state trajectories
%   x_v      - Noisy state trajectories
%   utraj    - Control input trajectories
%   dim_x    - State dimension
%   steps    - Number of time steps per trajectory
%   total_samples - Total number of columns
%
% Outputs:
%   U_data   - Input sequence matrix (U_-)
%   Y_0T     - Past output data (Y_-)
%   Y_1T     - Future output data (Y_+)

    x_meas_vec_0_v = zeros(dim_x, total_samples);
    x_meas_vec_1_v = zeros(dim_x, total_samples);
    x_meas_vec_0 = zeros(dim_x, total_samples);
    x_meas_vec_1 = zeros(dim_x, total_samples);
    u_mean_vec_0 = zeros(1, total_samples);

    idx_0 = 1;
    idx_1 = 1;

    for j = 1:dim_x:dim_x * total_samples / steps
        for i = 2:steps+1
            x_meas_vec_1_v(:,idx_1) = x_v(j:j+dim_x-1,i);
            x_meas_vec_1(:,idx_1)   = x(j:j+dim_x-1,i);
            idx_1 = idx_1 + 1;
        end

        for i = 1:steps
            u_mean_vec_0(:,idx_0)   = utraj(j,i);
            x_meas_vec_0(:,idx_0)   = x(j:j+dim_x-1,i);
            x_meas_vec_0_v(:,idx_0) = x_v(j:j+dim_x-1,i);
            idx_0 = idx_0 + 1;
        end
    end

    U_data = u_mean_vec_0(:,1:total_samples);
    Y_0T   = x_meas_vec_0_v(:,1:total_samples);
    Y_1T   = x_meas_vec_1_v(:,1:total_samples);
end
