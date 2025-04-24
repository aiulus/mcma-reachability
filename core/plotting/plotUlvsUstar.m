function plotUlvsUstar(time, ul_hist, u_hist, sys, output_dir, baseName)
% Plots learning inputs vs filtered (safe) inputs with input constraints
%
% Arguments:
%   time       - 1×T vector of time steps
%   ul_hist    - m×T matrix of learning (suggested) inputs
%   u_hist     - m×T matrix of safe (filtered) inputs
%   sys        - system struct with sys.dims.m and sys.bcs.U or sys.constraints.U
%   output_dir - output directory path for saving plots
%   baseName   - name prefix for saving files

    m = sys.dims.m;
    input_bounds = sys.bcs.U.interval;  % assume U has 'interval' property

    plots_per_fig = 3;
    num_figs = ceil(m / plots_per_fig);

    for f = 1:num_figs
        fig = figure('Position', [100, 100, 800, 600]);
        tiledlayout(plots_per_fig, 1);

        for sub = 1:plots_per_fig
            i = (f-1)*plots_per_fig + sub;
            if i > m, break; end

            nexttile;
            hold on; grid on;
            stairs(time, ul_hist(i,:), 'r:', 'LineWidth', 1.75, 'DisplayName', sprintf('$u_l^{(%d)}$', i));
            stairs(time, u_hist(i,:), 'b-', 'LineWidth', 1.25, 'DisplayName', sprintf('$u^{(%d)}$', i));

            % Bounds
            if ~isempty(input_bounds)
                lb = input_bounds(1,i);
                ub = input_bounds(2,i);
                yline(lb, 'm--', 'LineWidth', 1.25, 'DisplayName', 'Lower Bound');
                yline(ub, 'm--', 'LineWidth', 1.25, 'DisplayName', 'Upper Bound');
            end

            xlabel('$t$', 'Interpreter', 'latex');
            ylabel(sprintf('$u^{(%d)}$', i), 'Interpreter', 'latex');
            legend('Interpreter', 'latex');
            title(sprintf('Learning vs Safe Input %d', i), 'Interpreter', 'latex');
        end

        sgtitle('Learning vs Filtered Inputs', 'Interpreter', 'latex');

        % Save figure
        suffix = sprintf('%s_input_comp_fig%d', baseName, f);
        filename_pdf = fullfile(output_dir, [suffix, '.pdf']);
        filename_png = fullfile(output_dir, [suffix, '.png']);
        exportgraphics(fig, filename_pdf, 'ContentType','vector');
        exportgraphics(fig, filename_png, 'Resolution',600);
        close(fig);
    end
end


