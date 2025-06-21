function sys_approx = custom_conform_black(params,options,type)
% conform_black - identify a black-box model with genetic programming
%
% Syntax:
%    params = conform_black(params,options,type)
%
% Inputs:
%    params - parameters defining the conformance problem
%    options - options for the conformance checking
%    type - type of the algorithm
%
% Outputs:
%    params - parameters solving the conformance problem
%    fval - conformance cost
%    p_opt - estimated parameters
%    sys_upd - system object with the estimated parameters
%
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System
%        Identification," arXiv, 2024. 
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       31-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default parameters and options
sys = contDynamics();
%[params,options] = validateOptions(sys, params, options);

% reformat training and validation data
options.approx.p = 1;
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
sys_approx = nonlinearARX(type,f,dt, dim_y,dim_u, options.approx.p);
if options.approx.save_res
    save(options.approx.filename + "_sys", 'sys_approx');
end

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
% Corrected aux_testSuite2regress function
% x = [y_1(k-p) ... y_n(k-p) ... y_1(k-1) ... y_n(k-1) ...
%      u_1(k-p) ... u_n(k-p) ... u_1(k) ... u_n(k)]

if isempty(testSuite)
    x = [];
    y = [];
    return;
end

% Extract dimensions from the first test case
y_first = testSuite{1}.y;
u_first = testSuite{1}.u;
dim_y = size(y_first, 2);
dim_u = size(u_first, 2);
n_k_first = size(y_first,1);
n_s_first = size(y_first,3);

% Calculate expected number of features and total size
n_features = dim_y*p + dim_u*(p+1);
total_size = length(testSuite) * (n_k_first - p) * n_s_first;

x = zeros(total_size, n_features);
y = zeros(total_size, dim_y);
index = 1;

for m = 1:length(testSuite)
    y_m = testSuite{m}.y;
    u_m = testSuite{m}.u;
    n_s = size(y_m, 3); % number of samples in this test case

    for k = p+1:size(y_m,1)
        
        % Correctly reshape y_features per sample
        y_features = reshape(permute(y_m(k-p:k-1,:,:), [2 1 3]), 1, [], n_s);
        
        % Correctly reshape u_features per sample
        % This handles cases where u_m is 2D (same u for all samples) or 3D
        if ndims(u_m) == 2
            u_features_single = reshape(permute(u_m(k-p:k,:), [2 1]), 1, []);
            u_features = repmat(u_features_single, [1, 1, n_s]);
        else
            u_features = reshape(permute(u_m(k-p:k,:,:), [2 1 3]), 1, [], n_s);
        end
        
        % Concatenate features
        x_k = [y_features, u_features];

        y_k = y_m(k,:,:);

        % Assign to the main matrices
        end_index = index + n_s - 1;
        x(index:end_index, :) = squeeze(x_k)';
        y(index:end_index, :) = squeeze(y_k)';
        index = end_index + 1;
    end
end
end


% ------------------------------ END OF CODE ------------------------------