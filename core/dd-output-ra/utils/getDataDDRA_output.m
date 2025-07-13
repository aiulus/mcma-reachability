function [ytraj, ytraj_v, utraj] = getDataDDRA_output(sys_d, initpoints, steps, X0, W, V, u_seq)
    x0 = X0.c;
    ytraj(:, 1) = x0;
    index = 1;

    % System dimensions, assuming p = n (CORA)
    %dim_x = sys_d.nrOfStates;

    % System dimensions, assuming p = n (ss)
    dim_x = size(sys_d.A, 1);

    u = u_seq;

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