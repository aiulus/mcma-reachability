function [params, results] = custom_conform_black(params, options, type)
    % set default parameters and options
    sys_init = contDynamics();
    [params,options] = validateOptions(sys_init, params, options);
    
    % reformat training and validation data
    [xtrain,ytrain] = aux_testSuite2regress(params.testSuite_train, options.approx.p);
    [xval,yval] = aux_testSuite2regress(params.testSuite_val, options.approx.p);

    if options.approx.verbose
        fprintf("Approximate dynamics with %s. \n", type);
    end
    
    % run approximation
    switch type
        case "blackGP"
            f = aux_blackId_gp(xtrain, ytrain, xval, yval, params, options);
        case "blackCGP"
            f = aux_blackId_cgp(xtrain, ytrain, xval, yval, params, options);
    end

    % transform nonlinear function to sys object
    dim_y = size(params.testSuite{1}.y, 2);
    dim_u = size(params.testSuite{1}.u, 2);
    dt = params.testSuite{1}.sampleTime;
    sys = nonlinearARX(type,f,dt, dim_y,dim_u, options.approx.p);
    if options.approx.save_res
        save(options.approx.filename + "_sys", 'sys');
    end

    %% Preprocess testSuite for reachability
    params.testSuite = aux_preprocessTestSuite(sys, params.testSuite);
    params = aux_updateUWV(params); % Re-split U into U/W/V

    %% Compute reachset from learned model
    options_reach = options.options_reach;
    R = reach(sys, params.R0, params.U, options_reach);

    % Assign outputs    
    results.sys = sys;  
    params.sys = sys;
    results.R = R;
    results.fval = NaN;
    results.p = p_opt;
    results.unifiedOutputs = [];
end

% Auxiliary functions -----------------------------------------------------

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
    gp = rungp(@(x) config_gp(x, params, options, 'blackGP'));
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


function [x,y] = aux_testSuite2regress(testSuite, p)
    % x = [y_1(k-p) y_2(k-p) ... y_n(k-p) y_1(k-1) y_2(k-1) ... y_n(k-1) ...
    %      u_1(k-p) u_2(k-p) ... u_n(k-p) u_1(k) u_2(k) ... u_n(k)]
    
    total_size = length(testSuite) * (size(testSuite{1}.y,1)-p) * size(testSuite{1}.y, 3);
    x = zeros(total_size, size(testSuite{1}.y, 2)*p + size(testSuite{1}.u, 2)*(p+1));
    y = zeros(total_size, size(testSuite{1}.y, 2));
    index = 1;
    for m = 1:length(testSuite)
        y_m = testSuite{m}.y;
        u_m = testSuite{m}.u;
        for k = p+1:size(y_m,1)
            x_k = [reshape(permute(y_m(k-p:k-1,:,:), [2 1 3]), 1, [], size(y_m, 3)) ...
                repmat(reshape(permute(u_m(k-p:k,:), [2 1 3]), 1, [], 1), 1, 1, size(y_m, 3))];
            y_k = y_m(k,:,:);
            x(index:index+size(y_k, 3)-1, :) = squeeze(x_k)';
            y(index:index+size(y_k, 3)-1, :) = squeeze(y_k)';
            index = index + size(y_k, 3);
        end
    end
end
