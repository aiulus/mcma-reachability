function [U_data, Y_0T, Y_1T] = data2vec(u_traj, x_v, x, initpoints, n, steps)
    totalsamples = initpoints * steps;
    
    % prepeare Y_+ Y_-
    index_0 = 1;
    index_1 = 1;

    for j=1:n:initpoints*n
        for i=2:steps+1
            x_meas_vec_1_v(:,index_1) = x_v(j:j+n-1,i);
            x_meas_vec_1(:,index_1) = x(j:j+n-1,i);
            index_1 = index_1 +1;
        end
        for i=1:steps
            u_mean_vec_0(:,index_0) = u_traj(j,i);
            x_meas_vec_0(:,index_0) = x(j:j+n-1,i);
            x_meas_vec_0_v(:,index_0) = x_v(j:j+n-1,i);
            index_0 = index_0 +1;
        end
    end
    % U_data is U_-, Y_0T is Y_- , Y_1T is Y_+
    U_data = u_mean_vec_0(:,1:totalsamples); %same as u
    Y_0T = x_meas_vec_0_v(:,1:totalsamples);
    Y_1T = x_meas_vec_1_v(:,1:totalsamples);
end

