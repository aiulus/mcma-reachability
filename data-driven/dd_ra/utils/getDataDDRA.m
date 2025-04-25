function [x, utraj] = getDataDDRA(sys, initpoints, steps, X0, U, W)
    sys_d = sys.discrete;
    n = sys.dims.n; m = sys.dims.m;
    totalsamples = initpoints*steps;
    
    % Preallocate u
    u = zeros(m, totalsamples);
    for i=1:totalsamples
        u(:, i) = randPoint(U);
    end
    
    % Preallocate x, utraj
    x = zeros(initpoints*n, steps+1);
    utraj = zeros(initpoints*m, steps);

    % Initialize at X0.center
    x(1:n,1) = X0.center;
    
    index=1;
    for j=1:n:initpoints*n
        x(j:j+n-1,1) = randPoint(X0);
        for i=1:steps
            utraj(j:j+m-1,i) = u(:, index);
            x(j:j+n-1,i+1) = sys_d.A*x(j:j+n-1,i) + sys_d.B*u(index) + randPoint(W);      
            index=index+1;
        end
    end
end

