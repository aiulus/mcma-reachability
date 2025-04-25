function plotConstraints(time, y_hist, yl_hist, sys, output_dir, baseName)
% Plots outputs and filtered outputs with output bounds (e.g., from safety filter)
%
% Arguments:
%   time       - 1×T vector of time steps
%   y_hist     - p×T matrix of true system outputs
%   yl_hist    - p×T matrix of filtered outputs from the learning controller
%   sys        - struct with sys.dims.p and output bounds in sys.bcs.Y
%   output_dir - directory to save plots
%   baseName   - name prefix for output files

    p = sys.dims.p;
    output_bounds = sys.bcs.Y.interval;  % assume a 2×p matrix [inf; sup]

    plots_per_fig = 3;
    num_figs = ceil(p / plots_per_fig);

    for f = 1:num_figs
        fig = figure('Position', [100, 100, 800, 600]);
        tiledlayout(plots_per_fig, 1);

        for sub = 1:plots_per_fig
            i = (f-1)*plots_per_fig + sub;
            if i > p, break; end

            nexttile; hold on; grid on;

            % Plot y and yl
            stairs(time, y_hist(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('$y_{%d}$', i));
            stairs(time, yl_hist(i,:), 'r:', 'LineWidth', 1.5, 'DisplayName', sprintf('$y_l^{(%d)}$', i));

            % Output bounds
            lb = output_bounds(1,i);
            ub = output_bounds(2,i);
            yline(lb, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Lower Bound');
            yline(ub, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Upper Bound');

            xlabel('$t$', 'Interpreter', 'latex');
            ylabel(sprintf('$y_{%d}$', i), 'Interpreter', 'latex');
            legend('Interpreter', 'latex');
            title(sprintf('Output %d and Constraint Compliance', i), 'Interpreter', 'latex');
        end

        sgtitle('Output Trajectories with Constraints', 'Interpreter', 'latex');

        % Save
        suffix = sprintf('%s_output_comp_fig%d', baseName, f);
        filename_pdf = fullfile(output_dir, [suffix, '.pdf']);
        filename_png = fullfile(output_dir, [suffix, '.png']);
        exportgraphics(fig, filename_pdf, 'ContentType','vector');
        exportgraphics(fig, filename_png, 'Resolution',600);
        close(fig);
    end
end


