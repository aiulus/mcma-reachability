function f = aux_blackId_cgp(xtrain, ytrain, xval, yval, params, options)
    % genetic programming with conformance cost
    
    options.approx.xval = xval;
    options.approx.xtrain = xtrain;
    options.approx.yval = yval;
    options.approx.ytrain = ytrain;
    options.approx.pop_pre = true;
    tic
    % create initial propulation with normal genetic programming
    del_res = false;
    if ~isfield(options.approx, 'cgp_file_pop_pre') && ...
            options.approx.gp_num_gen > options.approx.cgp_num_gen
        options_gp = options;
        options_gp.approx.gp_num_gen = options.approx.gp_num_gen - options.approx.cgp_num_gen;
        options_gp.approx.gp_runs = 1;
        options_gp.approx.save_res = true;
        [~, fitness] = aux_blackId_gp(xtrain, ytrain, xval, yval, params, options_gp);
        options.approx.cgp_file_pop_pre = options.approx.filename;
        options.approx.cgp_conf_value = 5*fitness;
        if ~options.approx.save_res
            del_res = true;
        end
    end
    
    params.testSuite = params.testSuite_val;
    gp = rungp(@(x)config_gp(x, params, options, 'blackCGP'));
    T = toc;
    if options.approx.save_res
        save(options.approx.filename,"gp", "T");
    end
    if del_res
        % delete intermediary results
        for i_y = 1:size(ytrain,2)
            file_approx_iy = options.approx.filename + sprintf("_dim%d",i_y);
            delete(file_approx_iy+".mat");
        end
    end
    
    expr = gppretty(gp,'best');
    func = "@(y,u) [";
    for i_y = 1:size(ytrain,2)
        func = func + string(expr(i_y)) + ";";
    end
    f = eval(func + "]");
end