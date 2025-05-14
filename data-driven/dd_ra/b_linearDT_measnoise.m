% t_linear_measnoise - computes the data driven reachable set of discrete
% time systems with modeling and measurement noise
% x(k+1) = Ax(k) + Bu(k) + w(k)
% \tilde(x)(k) = x(k) + v(k)
%
% This example can be found in [1].
%
% Syntax:  
%    example_measnoise
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% References:
%    [1] Amr Alanwar, Anne Koch, Frank Allgöwer,Karl Henrik Johansson
%   "Data-Driven Reachability Analysis from Noisy Data"
% Author:       Amr Alanwar
% Written:      17-Jan-2021
% Last update:  21-May-2021 
% Last revision:---
%
% Adapted by:   Aybüke Ulusarslan
% Last update:  26-April-2025

%------------- BEGIN CODE --------------

rng(1);
clearvars;

%% Initial setup
% System initialization
systype = 'example0';
dt = 0.05;
sys = systemsDDRA(systype, dt);
sys_c = sys.cont;
sys_d = sys.discrete;
n = sys.dims.n;

% Catchall lookup dictionary
lookup = struct( ...
    'sys', sys, ...
    'initpoints', 10, ... % #trajectories
    'steps', 2, ... % #time steps
    'X0_center', ones(n, 1), ...
    'X0_spread', 0.1, ...
    'U_center', 10, ...
    'U_spread', 0.25, ...
    'W_center', 0, ...
    'W_spread', 0.005, ...
    'V_center', 0, ...
    'V_spread', 0.002 ...
    );

% Zonotope constructions
lookup.totalsamples = lookup.initpoints*lookup.steps;
lookup.X0 = zonotope(lookup.X0_center, lookup.X0_spread*diag(ones(n, 1))); % Initial set of states
lookup.U = zonotope(lookup.U_center, lookup.U_spread); % Set of inputs
lookup.W = zonotope(lookup.W_center*ones(n, 1), lookup.W_spread*ones(n, 1)); % Process noise bound
lookup.V = zonotope(lookup.V_center*ones(n, 1), lookup.V_spread*ones(n, 1)); % Measurement noise bound
lookup.Wmatzono= getGW(lookup);
lookup.Vmatzono = getGV(lookup);
lookup.AVmatzono = sys_d.A * lookup.Vmatzono;

%% Simulate & get trajectories
[U_full, X_0T, X_1T, X_0T_pure] = getTrajsMeasnoiseDDRA(lookup);
plot_toggle = 1;
[AB_av, AB_cmz] = estimateABmeasnoise(U_full, X_0T, X_1T, lookup, plot_toggle);

%% Plot trajectories
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;

% ------------------------ Refactor ------------------------
[AB, AV_oneterm] = estimate_AB_no_AV(X_1T, X_0T, U_full, Wmatzono, Vmatzono, W, V);

% Verify true AB inside set
AV_minus = check_true_AB_within(X_1T, AB, X_0T, U_full, Wmatzono, Vmatzono, sys_d, X_0T_pure);

% Compute reachable sets
totalsteps = 3;
redOrder = 390;
[X_model, X_data, X_data_av, X_data_cmz] = compute_next_sets(X0, U, W, V, sys_d, ...
    AB, AB_av, AB_cmz, AV_oneterm, totalsteps, redOrder);

% Visualization
projectedDims = {[1 2], [3 4], [4 5]};
axx = {[0.75, 1.5, 0.5, 4], [0.75, 3, 0.8, 2.2], [0.75, 2.3, 0.75, 2.8]};
visualize_reach_sets(X0, X_model, X_data, X_data_av, X_data_cmz, totalsteps, projectedDims, axx);


function [AB, AV_oneterm] = estimate_AB_no_AV(X_1T, X_0T, U_full, Wmatzono, Vmatzono, W, V)
    X1W_cen = X_1T - Wmatzono.center - Vmatzono.center;
    AB = X1W_cen * pinv([X_0T; U_full]);
    
    residual = X_1T - AB * ([X_0T; U_full]);
    minTerm = min(residual')';
    maxTerm = max(residual')';
    
    AV_oneterm = zonotope(interval(minTerm, maxTerm)) + (-1) * W + (-1) * V;
end

function AV_minus = check_true_AB_within(X_1T, AB, X_0T, U_full, Wmatzono, Vmatzono, sys_d, X_0T_pure)
    AV_minus_matzono = X_1T - AB * ([X_0T; U_full]) + (-1) * Wmatzono + (-1) * Vmatzono;
    AB_f = (AB * ([X_0T; U_full]) + AV_minus_matzono) * pinv([X_0T_pure; U_full]);
    
    intAB = intervalMatrix(AB_f).int;
    
    disp('True AB should be within estimated bounds:');
    disp('Upper bounds check:'), disp(all(intAB.sup >= [sys_d.A, sys_d.B], 'all'))
    disp('Lower bounds check:'), disp(all(intAB.inf <= [sys_d.A, sys_d.B], 'all'))
    
    VInt = intervalMatrix(AV_minus_matzono);
    leftLimit = VInt.Inf;
    rightLimit = VInt.Sup;
    AV_minus = zonotope(interval(min(leftLimit')', max(rightLimit')'));
end

function [X_model, X_data, X_data_av, X_data_cmz] = compute_next_sets(X0, U, W, V, sys_d, ...
    AB, AB_av, AB_cmz, AV_oneterm, totalsteps, redOrder)

    X_model = cell(totalsteps+1, 1); X_model{1} = X0;
    X_data = cell(totalsteps+1, 1); X_data{1} = X0;
    X_data_av = cell(totalsteps+1, 1); X_data_av{1} = X0;
    X_data_cmz = cell(totalsteps+1, 1); X_data_cmz{1} = conZonotope(X0);

    for i = 1:totalsteps
        % Model-based
        X_model{i+1} = sys_d.A * X_model{i} + sys_d.B * U + W;
        X_model{i+1} = reduce(X_model{i+1}, 'girard', redOrder);

        % Data-driven with AV
        X_data_av{i+1} = AB_av * (cartProd(X_data_av{i}, U)) + W;
        X_data_av{i+1} = reduce(X_data_av{i+1}, 'girard', redOrder);

        % Data-driven constrained matrix zonotope
        cart_cmz = cartProd(X_data_cmz{i}, conZonotope(U));
        X_data_cmz{i+1} = AB_cmz * cart_cmz + W;
        X_data_cmz{i+1} = reduce(X_data_cmz{i+1}, 'girard', redOrder);

        % Data-driven without AV
        if i == 1
            X_data{i+1} = AB * (cartProd(X_data{i}, U)) + AV_oneterm + W;
        else
            X_data{i+1} = AB * (cartProd(X_data{i} + V, U)) + AV_oneterm + W;
        end
        X_data{i+1} = reduce(X_data{i+1}, 'girard', redOrder);
    end
end

function visualize_reach_sets(X0, X_model, X_data, X_data_av, X_data_cmz, totalsteps, projectedDims, axx)

    numberofplots = totalsteps + 1;
    for plotRun = 1:length(projectedDims)
        figure('Renderer', 'painters', 'Position', [10 10 700 900]);
        hold on;
        
        % plot initial
        handleX0 = plot(X0, projectedDims{plotRun}, 'k-', 'LineWidth', 2);
        
        % plot reachable sets
        for iSet = 2:numberofplots
            handleModel = plot(X_model{iSet}, projectedDims{plotRun}, 'b', 'Filled', true, ...
                'FaceColor', [.8 .8 .8], 'EdgeColor', 'b');
        end
        for iSet = 2:numberofplots
            handleData = plot(X_data{iSet}, projectedDims{plotRun}, 'k-+');
        end
        for iSet = 2:numberofplots
            handleData_av = plot(X_data_av{iSet}, projectedDims{plotRun}, 'k');
        end
        for iSet = 2:numberofplots
            handleData_cmz = plot(X_data_cmz{iSet}, projectedDims{plotRun}, 'r', 'Template', 128);
        end
        
        xlabel(['x_', num2str(projectedDims{plotRun}(1))]);
        ylabel(['x_', num2str(projectedDims{plotRun}(2))]);
        
        legend([handleX0, handleModel, handleData, handleData_av, handleData_cmz], ...
            'Initial Set', 'Model-Based', 'Data-Driven', 'Data-Driven AV', 'Data-Driven CMZ', ...
            'Location', 'northwest');
        
        set(gca, 'FontSize', 20);
    end
end

