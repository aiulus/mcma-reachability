function [x, x_free, u] = getDataNonlinDDRA(sys, fun, initpoints, steps, X0, U, W)
    n = sys.dims.n; m = sys.dims.m;
    totalsamples = initpoints * steps;
    
    u = zeros(m, totalsamples);
    %input random sample points
    for i=1:totalsamples
        u(:,i) = randPointExtreme(params.U);
    end

    x = zeros(n * initpoints, steps + 1);
    x_free = zeros(n * initpoints, steps + 1);

    index = 1;
    for j = 1:n:initpoints * n
        x(j:j+n-1, 1) = randPoint(X0);
        x_free(j:j+n-1, 1) = x(j:j+n-1, 1);
        for i = 1:steps
            u(:, index) = randPoint(U);
            x_free(j:j+n-1, i+1) = fun(x(j:j+n-1, i), u(:, index));
            x(j:j+n-1, i+1) = x_free(j:j+n-1, i+1) + randPoint(W);
            index = index + 1;
        end
    end
end
