function [U_full, X_0T, X_1T] = ...
    getTrajsDDRA(sys, initpoints, T, x, utraj, plot_toggle)
    
    % For systemsDDRA objects
    %n = size(sys.A, 1); m = size(sys.B, 2);

    % For CORA loadDynamics() objects
    n = sys.nrOfStates; m = sys.nrOfInputs;
    totalsamples = initpoints*T;

    % Preallocate matrices
    u_mean_vec_0 = zeros(m, totalsamples);     
    x_meas_vec_0 = zeros(n, totalsamples);     
    x_meas_vec_1 = zeros(n, totalsamples);   

    index_0 = 1; index_1 = 1;

    for j = 1:n:initpoints * n
        for i = 2:T+1
            x_meas_vec_1(:, index_1) = x(j:j+n-1, i);
            index_1 = index_1 + 1;
        end
        for i = 1:T
            u_mean_vec_0(:, index_0) = utraj(j:j+m-1, i);
            x_meas_vec_0(:, index_0) = x(j:j+n-1, i);
            index_0 = index_0 + 1;
        end
    end
    
    % X_+ is X_1T
    % X_- is X_0T
    U_full = u_mean_vec_0(:,1:totalsamples); %same as u 
    X_0T = x_meas_vec_0(:,1:totalsamples);
    X_1T = x_meas_vec_1(:,1:totalsamples);
    
    %% Debug statement
    plot_toggle = 0;
    %%
    
    if plot_toggle
        figure;
        subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
        subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
        close;
    end
end

