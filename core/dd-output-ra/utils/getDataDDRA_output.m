function [ytraj, ytraj_v, utraj] = getDataDDRA_output(sys, initpoints, steps, X0, U, W)
    x0 = X0.c;
    ytraj(:, 1) = x0;
    index = 1;

    dim_x = sys.nrOfStates;

    for i=1:totalsamples
        u(i) = randPoint(U);
    end

    for j=1:dim_x:initpoints*dim_x
        ytraj(j:j+dim_x-1,1) = randPoint(X0); % initpoints many starting points sampled from X0
        ytraj_v(j:j+dim_x-1,1) =  ytraj(j:j+dim_x-1,1) + randPoint(V);
        
        for i=1:steps
            utraj(j,i) = u(index);
            ytraj(j:j+dim_x-1,i+1) = sys_d.A*ytraj(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);
            ytraj_v(j:j+dim_x-1,i+1) =  ytraj(j:j+dim_x-1,i+1) + randPoint(V);
            index=index+1;
        end
    end
end