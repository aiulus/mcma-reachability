clear; clc;
warning('off'); % Disabling warnings to avoid overwhelming the command 
                % window with CORA depreciation warnings --
                % TODO: Take care of the outdated usages causing them

%% Define experiment settings
systype = 'chain_of_integrators';
%dim_list = 2:2:10;
%noise_scales = [0.0, 0.1, 0.5, 1.0, 1.5, 2];
dim_list = [1, 5, 10];
noise_scales = [0.0, 1.0, 2];
plot_toggle = true;

results = struct();

for d = 1:length(noise_scales)
    dim = dim_list(d);

    for i = 1:length(noise_scales)
        scale_pn = noise_scales(i);
        for j = 1:length(noise_scales)
            scale_mn = noise_scales(j);
            noise_scale = struct('pcs', scale_pn, 'msmt', scale_mn);

            fprintf(['Running ZPC-DDSF for dim=%d, process noise x %.2f, ' ...
                'measurement noise x %.2f'], dim, scale_pn, scale_mn);
            tic;
            out = func_zonoDDSF(systype, dim, noise_scale, plot_toggle);
            total_time = toc;

            results.time(d, i, j) = total_time;
            results.cost(d, i, j) = out.Cost;
            results.cost_model(d, i, j) = out.Cost_model;
            results.conservatism(d, i, j) = out.conservatism;
        end
    end
end

%% Plot timing results
for i = 1:length(noise_scales)
    figure;
    for j = 1:length(noise_scales)
        plot(dim_list, results.time(:, i, j), '-o', 'DisplayName',sprintf('Mes. Noise %.2f', noise_scales(j)));
        hold on;
    end
    xlabel('System Dimension'); ylabel('Computation Time [s]');
    title(sprintf('ZPC-DDSF Runtime (Process Noise = %.2f)', noise_scales(i)));
    legend; grid on;
end

%% Plot control cost
for i = 1:length(noise_scales)
    figure;
    for j = 1:length(noise_scales)
        plot(dim_list, results.cost(:,i,j), '-o', 'DisplayName', sprintf('ZPC Cost, Mes. %.2f', noise_scales(j)));
        hold on;
        plot(dim_list, results.cost_model(:,i,j), '--s', 'DisplayName', sprintf('RMPC Cost, Mes. %.2f', noise_scales(j)));
    end
    xlabel('System Dimension'); ylabel('Control Cost');
    title(sprintf('ZPC vs RMPC Cost (Meas Noise = %.2f)', noise_scales(i)));
    legend; grid on;
end

%% Plot conservatism
for i = 1:length(noise_scales)
    figure;
    for j = 1:length(noise_scales)
        plot(dim_list, results.conservatism(:, i, i), '-*', 'DisplayName', sprintf('Mes. %.2f',noise_scales(j)));
        hold on;
    end
    xlabel('System Dimension'); ylabel('Conservatism');
    title(sprintf('ZPC Conservatism (Proc. Noise = %.2f)', noise_scales(i)));
    legend; grid on;
end
warning('on');