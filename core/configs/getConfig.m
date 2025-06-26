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
