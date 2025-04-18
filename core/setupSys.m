function [sys_d, ref, intc, y0, U, X0, W, V, VAmatzono, samplingtime] = setupSys()
    dim_x = 5;
    A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
    B_ss = ones(5,1);
    C = [1,0,0,0,0]; D = 0;
    sys_c = ss(A, B_ss, C, D);

    samplingtime = 0.05;
    sys_d = c2d(sys_c, samplingtime);

    uref = 8;
    ref = inv(eye(5)-sys_d.A) * sys_d.B * uref;

    y_lb = [-10; 2; -10; -10; -10]; 
    y_ub = [10; 10; 10; 10; 10]; 
    intc = interval(y_lb, y_ub);

    y0 = [-2; 4; 3; -2.5; 5.5];
    X0 = zonotope([y0, 25 * diag(ones(5,1))]);
    U = zonotope([uref - 1, 20 - 1]);

    wfac = 0.01; 
    W = zonotope([zeros(5,1), wfac * ones(5,1)]);

    vfac = 0.002; 
    V = zonotope([zeros(5,1), vfac * ones(5,1)]);

    VAmatzono = sys_d.A * V;
end
