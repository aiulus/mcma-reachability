function [x_meas_vec_0, x_meas_vec_1, x_free_vec_0, x_free_vec_1, ...
    U_full, X_0T, X_1T] = getTrajsNonlinDDRA(lookup)

    U = lookup.U; R0 = lookup.R0; W = lookup.W;  
    n = lookup.n; m = lookup.sys.dims.m; 
    fun = lookup.fun; 
    initpoints = lookup.initpoints; 
    steps = lookup.steps;
    totalsamples = lookup.totalsamples;

    % Preallocate
    u = zeros(m, totalsamples);
    x = zeros(n * initpoints, steps + 1);
    x_free = zeros(n * initpoints, steps + 1);
    x_meas_vec_0 = zeros(n, totalsamples);
    x_meas_vec_1 = zeros(n, totalsamples);
    x_free_vec_0 = zeros(n, totalsamples);
    x_free_vec_1 = zeros(n, totalsamples);
    
    % Sample random input vectors from U
    for i=1:totalsamples
        % u(:,i) = randPointExtreme(U);
        u(:,i) = randPoint(U);
    end
    
    % Simulate trajectories
    index = 1;
    for j=1:n:initpoints*n
        x(j:j+n-1,1) = randPoint(R0);
        x_free(j:j+n-1,1) = x(j:j+n-1,1);
        for i=1:steps
            x_free(j:j+n-1,i+1) = fun(x(j:j+n-1,i),u(:,index));
            x(j:j+n-1,i+1) = fun(x(j:j+n-1,i),u(:,index)) +randPoint(W);
            index = index + 1;
        end
    end    
    
    % Combine trajectories
    index_0 =1;
    index_1 =1;
    for j=1:n:initpoints*n
        for i=2:steps+1        
            x_meas_vec_1(:,index_1) = x(j:j+n-1,i);
            x_free_vec_1(:,index_1) = x_free(j:j+n-1,i);
            index_1 = index_1 +1;
        end
        for i=1:steps
            x_free_vec_0(:,index_0) = x_free(j:j+n-1,i);
            x_meas_vec_0(:,index_0) = x(j:j+n-1,i);
            index_0 = index_0 +1;
        end
    end
    U_full = u(:,1:totalsamples);
    X_0T = x_meas_vec_0(:,1:totalsamples);
    X_1T = x_meas_vec_1(:,1:totalsamples);
end
