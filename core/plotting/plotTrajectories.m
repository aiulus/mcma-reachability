function plotTrajectories(Rplotall, y_t, YPred, dims, output_dir, baseName)
% Plots reachable sets R_t along with the actual system trajectory y_t
% and predicted nominal trajectory YPred.
%
% Arguments:
%   Rplotall   - cell array of zonotopes or intervals (reach sets at each t)
%   y_t        - noisy system output, size [n x T+1]
%   YPred      - predicted output, size [n x T+1]
%   dims       - vector of dimensions to plot, e.g., [1 2] or [1 2 3]
%   output_dir - directory where plots should be saved
%   baseName   - name prefix for file naming

    assert(length(dims) == 2 || length(dims) == 3, 'dims must be 2D or 3D');

    maxsteps = length(Rplotall);

    fig = figure('Renderer', 'painters', 'Position', [100, 100, 900, 800]);
    hold on;
    box on;

    % Plot reachable sets
    for k = 1:maxsteps
        if ~isempty(Rplotall{k})
            plot(Rplotall{k}, dims, 'r', 'Filled', false);
        end
    end

    % Plot actual noisy trajectory
    if length(dims) == 3
        plot3(y_t(dims(1), 2:end), y_t(dims(2), 2:end), y_t(dims(3), 2:end), '+b');
        plot3(YPred(dims(1), 2:end), YPred(dims(2), 2:end), YPred(dims(3), 2:end), '*k');
    else
        plot(y_t(dims(1), 2:end), y_t(dims(2), 2:end), '+b');
        plot(YPred(dims(1), 2:end), YPred(dims(2), 2:end), '*k');
    end

    legend({'Reachable set $\hat{\mathcal{R}}_t$', ...
            'System trajectory $y(t)$', ...
            'Predicted trajectory $\hat{y}(t)$'}, ...
            'Interpreter','latex','Location','best');
        
    xlabel(['$y_{', num2str(dims(1)), '}(t)$'], 'Interpreter', 'latex');
    ylabel(['$y_{', num2str(dims(2)), '}(t)$'], 'Interpreter', 'latex');
    if length(dims) == 3
        zlabel(['$y_{', num2str(dims(3)), '}(t)$'], 'Interpreter', 'latex');
    end

    set(gca, 'FontSize', 16);
    grid on;
    
    % Save plot
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    suffix = sprintf('%s_dims_%s', baseName, strjoin(string(dims), '_'));
    filename_pdf = fullfile(output_dir, [suffix, '.pdf']);
    filename_png = fullfile(output_dir, [suffix, '.png']);
    
    exportgraphics(fig, filename_pdf, 'ContentType','vector');
    exportgraphics(fig, filename_png, 'Resolution',600);
    close(fig);
end

