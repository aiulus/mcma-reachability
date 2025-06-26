function out = func_conformance_linear(systype, dim, noise_scale, plot_toggle)
    rng(1);

    dt = 0.05;
    sys_full = systemsDDRA(systype, dt, dim);
    sys = sys_full.CORA;

    X0_center = ones(sys.nrOfStates, 1);
    X0_spread = 0.1;
    U_center = 10;
    U_spread = 0.25;
    W_spread = 0.005 * noise_scale;

    R0 = zonotope([X0_center, X0_spread * eye(sys.nrOfStates)]);
    U_set = zonotope([U_center * ones(sys.nrOfInputs, 1), ...
                     U_spread * eye(sys.nrOfInputs)]);

    %% Build params
    params = struct;
    params.R0 = R0;
    params.U = U_set;
    params.tFinal = dt * 20;
    
    %% Create test suite
    options_test.p_extr = 0.2;
    options_test.inputCurve = "rand";
    params.testSuite = createTestSuite(sys, params, 20, 30, 1, options_test);

    %% Conformance options
    options.cs.robustnessMargin = 1e-9;
    options.cs.cost = 'interval';
    options.cs.verbose = false;

    %% Synthesize conformance
    tic;
    [params_conform, ~] = conform(sys, params, options);
    total_time = toc;

    %% Validate
    configs = {struct('sys', sys, 'params', params_conform, 'options', options, 'name', 'cora')};
    testSuite{1} = params.testSuite{1};
    check_contain = false;
    plot_settings = struct('dims', [1 2], 'name', 'CORA Eval', 's_val', 1);

    if plot_toggle
        validateReach(testSuite{1}, configs, check_contain, plot_settings);
    end

    %% Output
    out = struct;
    out.time = total_time;
    out.configs = configs;
    out.testSuite = testSuite;
end



