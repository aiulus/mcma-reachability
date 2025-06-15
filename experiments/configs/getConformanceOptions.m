function options = getConformanceOptions(options_reach, cost_norm, constraints, sys)
    options = options_reach;

    options.cs.robustnessMargin = 1e-9;
    options.cs.verbose = false;
    options.cs.cost = cost_norm;
    options.cs.constraints = constraints;
    
    % Black-box approximation options
    options.approx.gp_parallel = true;
    options.approx.gp_pop_size = 50;
    options.approx.gp_num_gen = 30;
    options.approx.gp_func_names = {'times','plus', 'square'};
    options.approx.gp_max_genes = 2;
    options.approx.gp_max_depth = 2;
    options.approx.gp_parallel = false;
    options.approx.cgp_num_gen = 5;
    options.approx.cgp_pop_size_base = 5;
    options.approx.save_res = false;
    %options.approx.p = sys.n_p;
    %options.approx.p = sys.nrOfOutputs;
    options.approx.p = 0;
end