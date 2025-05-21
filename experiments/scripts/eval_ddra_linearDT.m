clear; clc;

%% Define experiment settings
system_type = 'chain_of_integrators';
dim_list = 2:2:10; % System dimensions
noise_scales = [0.0, 0.1, 0.5, 1.0, 1.5, 2]; % Noise scaling factors
measure_type = 'volume'; % Or: 'bbox_volume', 'fro_norm'
plot_toggle = true;

results = struct();

%% Run adapted <a_linearDT.m> for each dimension and noise level
for d = 1:length(dim_list)
    dim = dim_list(d);

    for ns = 1:length(noise_scales)
        scale_noise = noise_scales(ns);

        % Run main function
        fprintf('Running DDRA for dim = %d, noise scale = %.2f\n', dim, scale_noise);
        tic;
        out = func_ddra_linearDT(system_type, dim, scale_noise, plot_toggle);
        total_time = toc;

        % Compute size metrics at the final time step
        lastSet = out.X_data{end};
        size_val = getReachsetSize(lastSet, measure_type);

        results.time(d, ns) = total_time;
        results.size(d, ns) = size_val;
    end
end

%% Plot the results
figure;
for ns = 1:length(noise_scales)
    plot(dim_list, results.time(:, ns), '-o', 'DisplayName', sprintf('Noise %.2f', noise_scales(ns)));
    hold on;
end
xlabel('System Dimension'); ylabel('Computation Time [s]');
title('<a_linearDT> Runtime vs. Dimension'); legend; grid on;

figure;
for ns = 1:length(noise_scales)
    plot(dim_list, results.size(:, ns), '-s', 'DisplayName',sprintf('Noise %s', noise_scales(ns)));
    hold on;
end
xlabel('System Dimension'); ylabel(['Reachset Size (' measure_type ')']);
title('Reachset Size vs. Dimension'); legend; grid on;
