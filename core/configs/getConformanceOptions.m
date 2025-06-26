function options = getConformanceOptions(options_reach, cost_norm, constraints, sys)
    options = options_reach;

    options.cs.robustnessMargin = 1e-9;
    options.cs.verbose = false;
    options.cs.cost = cost_norm;
    options.cs.constraints = constraints;
    options.cs.derivRecomputation = false;
    %% TODO: Remove / Debug statement
    %options.cs.w = zeros(100, 1);
    %options.cs.w = ones(1, settings.n_k);
    options.cs.w = ones(1, sys.nrOfStates);
    
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
    options.approx.p = 1;
    options.approx.verbose = false;
    options.approx.filename = 'temp_gp_results';
    options.approx.gp_runs = 10;
    %% Achtung! Not sure what this does
    options.approx.cgp_n_m_conf = 4;
end