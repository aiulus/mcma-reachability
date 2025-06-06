function completed = flexBlackBoxConform(varargin)
% --------------------------------------------------------------------
    % HOW TO USE flexBlackBoxConform:
    %
    %   flexBlackBoxConform()
    %     – Uses all default settings:
    %         • dynamics           = "Square"
    %         • plot_settings      = config.plot_settings from getConfig()
    %         • cost_norm          = "interval"
    %         • constraints        = "half"
    %         • identification     = ["blackGP","blackCGP"] + "true"
    %         • reachability opts  = config.options_reach from getConfig()
    %         • test‐suite params  = config.options_testS.p_extr = 0.3
    %
    %   flexBlackBoxConform('dynamics', DYN_NAME)
    %     – Override the system dynamics. DYN_NAME must be a valid input
    %       to loadDynamics (e.g. "VanDerPol", "Square", etc.).
    %
    %   flexBlackBoxConform('plot_settings', PLOT_STRUCT)
    %     – Override the default plotting settings. PLOT_STRUCT should have
    %       the same fields as returned by getConfig().plot_settings, namely:
    %         .plot_Yp  (logical, default=false)
    %         .dims     (1×2 vector, default=[1 2])
    %       If you omit this name-value pair, plot_settings is taken from
    %       getConfig().plot_settings.
    %
    %   flexBlackBoxConform('dynamics', DYN_NAME, 'plot_settings', PLOT_STRUCT)
    %     – Override both at once. Name-value pairs can appear in any order.
    %
    % IMPORTANT DEFAULTS (defined in getConfig()):
    %   • settings.n_m         = 2
    %   • settings.n_s         = 50
    %   • settings.n_k         = 4
    %   • settings.n_m_train   = 100
    %   • settings.n_s_train   = 10
    %   • settings.n_k_train   = 4
    %   • settings.n_m_val     = 5
    %   • settings.n_s_val     = 10
    %   • settings.n_k_val     = 4
    %   • options_reach.zonotopeOrder     = 100
    %   • options_reach.tensorOrder       = 2
    %   • options_reach.errorOrder        = 1
    %   • options_reach.tensorOrderOutput = 2
    %   • options_reach.verbose           = false
    %   • config.options_testS.p_extr      = 0.3
    %   • config.plot_settings.plot_Yp     = false
    %   • config.plot_settings.dims        = [1 2]
    %   • cost_norm  = "interval"
    %   • constraints = "half"
    %
    % To change any of the “non‐varargin” defaults (e.g. cost_norm, constraints,
    % or the numbers of trajectories), edit getConfig() or the hard-coded lines
    % inside this function directly.
    % --------------------------------------------------------------------
    rng(2)
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'plot_settings', []);
    addParameter(p, 'dynamics', "Square");
    parse(p, varargin{:});

    override_plot_settings = p.Results.plot_settings;
    dynamics = p.Results.dynamics;

    %% Load default config
    config = getConfig();
    settings = config.settings;
    options_reach = config.options_reach;
    options_testS.p_extr = 0.3;

    % Use default unless override is provided
    if isempty(override_plot_settings)
        plot_settings = config.plot_settings;
    else
        plot_settings = override_plot_settings;
    end
    
    cost_norm = "interval"; % norm for the reachable set: "interval","frob"    
    
    constraints = "half"; % constraint type: "half", "gen"    
    
    methodsGray = ["blackGP","blackCGP"]; % identification approach    
    methods = ["true" methodsGray];

    % Load system dynamics
    [sys, params_true.R0, params_true.U, p_true] = loadDynamics(dynamics, "rand");
    params_true.tFinal = sys.dt * settings.n_k - sys.dt;

    % Create identification data
    params_true.testSuite = createTestSuite(sys, params_true, settings.n_k, settings.n_m, settings.n_s, options_testS);
    params_true.testSuite_train = createTestSuite(sys, params_true, settings.n_k_train, settings.n_m_train, settings.n_s_train);
    params_true.testSuite_val = createTestSuite(sys, params_true, settings.n_k_val, settings.n_m_val, settings.n_s_val);

    %% Conformance Identification ---------------------------------------------
    % Get default identification and black-box approximation options
    options = getConformanceOptions(options_reach, cost_norm, constraints, sys);

    % Create struct for saving the identification results for each system
    %% TODO: rename to 'results'
    results = cell(length(methodsGray)+1, 1);
    results{1}.sys = sys;
    results{1}.params = rmfield(params_true, 'testSuite');
    results{1}.options = options_reach;
    results{1}.name = "true";

    % Initial Estimates of the Disturbance Sets
    c_R0 = center(params_true.R0);
    c_U = center(params_true.U);

    params_id_init = params_true;
    params_id_init.R0 = zonotope(c_R0);
    params_id_init.U = zonotope([c_U eye(size(c_U, 1)) ones(size(c_U))]);

    % Identification ------------------------------------------------------
    for i = 1:length(methodsGray)
        type = methodsGray(i);
        fprintf("Identification with method %s \n", type);

        tic;
        [results{i+1}.params, results] = conform(sys, params_id_init, options, type);
        Ts = toc;

        results{i+1} = struct( ...
            'sys', results.sys, ...
            'options', options_reach, ...
            'name', type ...
            );
        
        fprintf("Identification time: %.4f\n", Ts);
    end

    %% Validation and Visualization -------------------------------------------

    % Sanity check: Compute Reachable Sets and Check Containment of the 
    % Identification Test Cases    
    validateReachableSets(params_true.testSuite, results, settings.n_k_val, ...
        ["true" methodsGray], 'label', 'IDENTIFICATION DATA', 'require_full_containment', true);


    % Create Validation Data
    params_true.tFinal = sys.dt * settings.n_k_val - sys.dt;
    testSuite_val = createTestSuite(sys, params_true, settings.n_k_val, settings.n_m_val, ...
        settings.n_s_val, options_testS);

    % Compute Reachable Sets and Check Containment of the Validation Test Cases
    validateReachableSets(testSuite_val, results, settings.n_k_val, ...
        methods, 'plot_settings', plot_settings, 'label', 'VALIDATION DATA');
    
    % example completed
    completed = true;
end

function validateReachableSets(testSuite, configs, n_k_val, methods, varargin)
    % Parse optional inputs
    p = inputParser;
    addOptional(p, 'plot_settings', []);
    addOptional(p, 'label', 'VALIDATION DATA');
    addOptional(p, 'require_full_containment', false);
    parse(p, varargin{:});
    plot_settings = p.Results.plot_settings;
    label = p.Results.label;
    require_full_containment = p.Results.require_full_containment;

    num_out = 0;
    num_in = 0;
    check_contain = 1;

    for m = 1:length(testSuite)
        if isempty(plot_settings)
            [~, eval] = validateReach(testSuite{m}, configs, check_contain);
        else
            [~, eval] = validateReach(testSuite{m}, configs, check_contain, plot_settings);
        end
        num_out = num_out + eval.num_out;
        num_in = num_in + eval.num_in;
    end

    num_all = length(testSuite) * n_k_val * size(testSuite{1}.y, 3);
    fprintf("%s: \n", label);
    for i = 1:length(configs)
        p_contained = 100 - (num_out(i)/(num_out(i)+num_in(i))) * 100;
        suffix = require_full_containment * " (must be 100%!)";
        fprintf("%s: %.2f%% of the samples are contained in the reachable set%s. \n", ...
            methods(i), p_contained, suffix);
        fprintf("%s: %.2f%% of the samples were not valid. \n", ...
            methods(i), (num_all - (num_out(i) + num_in(i))) / num_all * 100);
    end
end

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
    options.approx.p = sys.n_p;
end

function config = getConfig()
    settings.n_m = 2; % #different input trajectories
    settings.n_s = 50; % #sample trajectories per input trajectory
    settings.n_k = 4; % Length of the identification trajectories

    % Training and validation data
    settings.n_m_train = 100;
    settings.n_s_train = 10;
    settings.n_k_train = 4;
    settings.n_m_val = 5;
    settings.n_s_val = 10;
    settings.n_k_val = 4;    

    config.settings = settings;

    % Reachability settings
    options_reach.zonotopeOrder = 100;
    options_reach.tensorOrder = 2;
    options_reach.errorOrder = 1;
    options_reach.tensorOrderOutput = 2;
    options_reach.verbose = false;

    config.options_reach = options_reach;

    config.options_testS.p_extr = 0.3;
    
    % Plotting
    config.plot_settings.plot_Yp = false;
    config.plot_settings.dims = [1 2];
end



