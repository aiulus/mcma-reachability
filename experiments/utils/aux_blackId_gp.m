function [f, fitness_sum] = aux_blackId_gp(xtrain, ytrain, xval, yval, params, options)
    % normal genetic programming
    
    options.approx.xval = xval;
    options.approx.xtrain = xtrain;
    fitness_sum = 0;
    func = "@(y,u) [";
    for i_y = 1:size(ytrain,2)
        options.approx.ytrain = ytrain(:,i_y);
        options.approx.yval = yval(:,i_y);
        file_approx_iy = options.approx.filename + sprintf("_dim%d",i_y);
    
        % start new gp run
        if options.approx.verbose
            fprintf('Output dimension %d. \n', i_y)
        end
        tic
        gp = rungp(@(x) custom_config_gp(x, params, options, 'blackGP'));
        T = toc;
        if options.approx.save_res
            save(file_approx_iy,"gp", "T");
        end
        expr = gppretty(gp,'valbest');
        func = func + string(expr) + ";";
        fitness_sum = fitness_sum + gp.results.best.fitness;
    end
    f = eval(func + "]");
end