function plotTerminalSet(RN, S_f, dims, output_dir, baseName)
% Plots the terminal reachable set R_N and the terminal safe set S_f
% in a specified 2D projection.
%
% Arguments:
%   RN         - zonotope R{N+1}, final reachable set
%   S_f        - zonotope terminal safe set
%   dims       - dimensions to project (e.g., [1 2])
%   output_dir - directory for saving figures
%   baseName   - string used as part of filename

    assert(length(dims) == 2, 'Only 2D projections are supported for terminal inclusion plot');

    fig = figure('Position', [100, 100, 700, 600]);
    hold on;
    box on;

    % Plot terminal safe set (background)
    plot(S_f, dims, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'k', 'LineStyle', '--');
    
    % Plot terminal reachable set R_N
    plot(RN, dims, 'r', 'LineWidth', 1.5);

    legend({'Terminal Safe Set $\mathcal{S}_f$', 'Reachable Set $R_N$'}, ...
           'Interpreter', 'latex', 'Location', 'best');

    xlabel(['$y_{', num2str(dims(1)), '}$'], 'Interpreter', 'latex');
    ylabel(['$y_{', num2str(dims(2)), '}$'], 'Interpreter', 'latex');
    title('Terminal Set Inclusion: $R_N \subseteq \mathcal{S}_f$', 'Interpreter', 'latex');

    set(gca, 'FontSize', 16);
    grid on;

    % Save figure
    suffix = sprintf('%s_terminal_inclusion_%s', baseName, strjoin(string(dims), '_'));
    filename_pdf = fullfile(output_dir, [suffix, '.pdf']);
    filename_png = fullfile(output_dir, [suffix, '.png']);
    exportgraphics(fig, filename_pdf, 'ContentType','vector');
    exportgraphics(fig, filename_png, 'Resolution',600);
    close(fig);
end
