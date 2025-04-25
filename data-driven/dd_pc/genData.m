function [u_traj, x_v, x] = genData(sys, X0, U, W, V, initpoints, totalsamples)
    sys_d = sys.discrete;
    n = size(sys_d.A, 1);
    steps = totalsamples / initpoints;

    % randomly choose constant inputs for each step / sampling time
    for i=1:totalsamples
        u(i) = randPoint(U);
    end
    
    % generate data from different trajectories with noise
    x0 = X0.center;
    x(:,1) = x0;

    index=1;
    for j=1:n:initpoints*n
        x(j:j+n-1,1) = randPoint(X0);
        x_v(j:j+n-1,1) =  x(j:j+n-1,1) + randPoint(V);
        
        for i=1:steps
            u_traj(j,i) = u(index);
            x(j:j+n-1,i+1) = sys_d.A*x(j:j+n-1,i) + sys_d.B*u(index) + randPoint(W);
            x_v(j:j+n-1,i+1) =  x(j:j+n-1,i+1) + randPoint(V);
            index=index+1;
        end
    end
end

