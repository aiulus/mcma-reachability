function ts_out = reshapeTestSuiteLinear(ts_in, sys)
    n_m = length(ts_in);
    ts_out = cell(n_m,1);
    sysClass = class(sys);
    
    for i=1:n_m
        ts_in_i = ts_in{i};
        y_i = ts_in_i.y;
        x0_i = ts_in_i.initialState;
        u_i = ts_in_i.u;
        ts_out{i} = testCase(y_i, u_i, x0_i, sys.dt, sysClass);
    end
end