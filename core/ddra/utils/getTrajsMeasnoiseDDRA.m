function [U_full, X_0T, X_1T, X_0T_pure] = getTrajsMeasnoiseDDRA(lookup, plot_toggle)
    %% Extract relevant parameters
    X0 = lookup.X0; U = lookup.U; W = lookup.W; V = lookup.V;
    n = lookup.sys.dims.n; m = lookup.sys.dims.m;
    initpoints = lookup.initpoints; steps = lookup.steps; totalsamples = lookup.totalsamples;
    sys_d = lookup.sys.discrete;

    %% Preallocation
    u = zeros(m, totalsamples);            % inputs
    x = zeros(n*initpoints, steps+1);       % nominal states
    x_v = zeros(n*initpoints, steps+1);     % measured states
    utraj = zeros(n*initpoints, steps);     % input trajectories

    x_meas_vec_1_v = zeros(n, totalsamples); % measured x(k+1)
    x_meas_vec_0_v = zeros(n, totalsamples); % measured x(k)
    x_meas_vec_0   = zeros(n, totalsamples); % true x(k)
    u_mean_vec_0   = zeros(m, totalsamples); % inputs u(k)


    for i=1:totalsamples
        u(i) = randPoint(U);
    end
    
    % simulate the discrete system starting from random initial points
    %x(:,1) = randPoint(X0);
    index=1;
    for j=1:n:initpoints*n
        x(j:j+n-1,1) = randPoint(X0);
        x_v(j:j+n-1,1) = x(j:j+n-1,1);
        for i=1:steps
            utraj(j,i) = u(index);
            x(j:j+n-1,i+1) = sys_d.A*x(j:j+n-1,i) + sys_d.B*u(index) + randPoint(W);  
            x_v(j:j+n-1,i+1) =  x(j:j+n-1,i+1) + randPoint(V);
            index=index+1;
        end
    end   
    
    index_0 =1;
    index_1 =1;
    for j=1:n:initpoints*n
        for i=2:steps+1
            x_meas_vec_1_v(:,index_1) = x_v(j:j+n-1,i);
            index_1 = index_1 +1;
        end
        for i=1:steps
            u_mean_vec_0(:,index_0) = utraj(j,i);
            x_meas_vec_0_v(:,index_0) = x_v(j:j+n-1,i);
            x_meas_vec_0(:,index_0) = x(j:j+n-1,i);
            index_0 = index_0 +1;
        end
    end

    U_full = u_mean_vec_0(:,1:totalsamples); 
    X_0T = x_meas_vec_0_v(:,1:totalsamples);
    X_1T = x_meas_vec_1_v(:,1:totalsamples);
    X_0T_pure = x_meas_vec_0(:,1:totalsamples);

    if plot_toggle == 1
        %% Plot trajectories
        figure;
        subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
        subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
        close;
    end
end

