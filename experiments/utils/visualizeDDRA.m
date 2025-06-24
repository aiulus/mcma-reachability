function visualizeDDRA(X0, X_model, X_data, projectedDims, axx, n_k, varargin)
% visualizeDDRA  –  2-D projection of DDRA reachable sets (model ⬄ data)
%
%   This is a **drop-in replacement** for the ad-hoc function
%   `visualizeAlinearDT`.  Internally it copies the plotting style that is
%   used by CORA’s `testCase/validateReach`, so figures look identical to
%   the white-/grey-/black-box pipeline produced by `validateReach`.
%
%   INPUTS
%   ======
%     X0            – initial zonotope (plotted once)
%     X_model{1..k} – cell array returned by *propagateDDRA* (model set)
%     X_data {1..k} – cell array returned by *propagateDDRA* (data set)
%     projectedDims – cell array, each entry contains a 2-vector with the
%                     output indices you want to plot, e.g. {[1 2],[3 4]}
%     axx           – (optional) axis limits per projection, same length
%                     as *projectedDims*; leave empty to auto-scale
%     n_k           – number of timesteps you want to display (usually
%                     `length(X_model)`)
%
%   OPTIONAL name–value pairs (same keys as validateReach)
%   -----------------------------------------------------
%     'colors'      – n×3 RGB matrix                (default = lines())
%     'lines'       – string array of line styles   (default = see code)
%     'linewidths'  – numeric vector                (default = 0.1)
%     'k_plot'      – timesteps to show             (default = 1:n_k)
%     'name'        – title prefix                  (default auto)
%
%   EXAMPLE
%   =======
%       projectedDims = {[1 2]};
%       visualizeDDRA(R0, X_model_P2, X_data_P2, projectedDims, {}, ...
%                     length(X_model_P2));
%
% -------------------------------------------------------------------------
%   © 2025  Your Name / Lab
% -------------------------------------------------------------------------

% ---------------- Input parsing ---------------------------------------- %

default.lines      = ["-", "--"];
default.colors     = lines(2);
default.linewidths = 0.1*ones(2,1);

a = inputParser;
a.KeepUnmatched = true;     % ignore unknown keys silently
addParameter(a,'colors',     default.colors);
addParameter(a,'lines',      default.lines);
addParameter(a,'linewidths', default.linewidths);
addParameter(a,'k_plot',     1:n_k);
addParameter(a,'name',       "DDRA reachability comparison");
parse(a, varargin{:});
ps = a.Results;                           % plot_settings clone

if isempty(projectedDims)
    projectedDims = {[1 2]};              % fallback
end
num_row = ceil( sqrt( numel(ps.k_plot) ) );
num_col = ceil( numel(ps.k_plot) / num_row );

% ---------------- Figure(s) -------------------------------------------- %
for i_dim = 1:numel(projectedDims)
    dim = projectedDims{i_dim};   % 2-vector
    fig = figure('Units','normalized','OuterPosition',[0 0 1 1], ...
                 'Name',sprintf('visualizeDDRA – dims [%d %d]',dim));

    for i_k = 1:numel(ps.k_plot)
        k      = ps.k_plot(i_k);
        subplot(num_col,num_row,i_k); hold on;

        % initial set (only once)
        if k == 1
            plot(X0, dim, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','k', ...
                 'LineWidth',2, 'DisplayName','$X_0$');
        end

        % model-based reachable set
        if k <= numel(X_model) && isa(X_model{k},'contSet')
            plot(X_model{k}, dim, 'FaceColor',ps.colors(1,:), ...
                 'FaceAlpha',0.1, 'EdgeColor',ps.colors(1,:), ...
                 'LineStyle',ps.lines(1), 'LineWidth',ps.linewidths(1), ...
                 'DisplayName','$\mathcal X_{model}$');
        end

        % data-driven reachable set
        if k <= numel(X_data) && isa(X_data{k},'contSet')
            plot(X_data{k}, dim, 'FaceColor',ps.colors(2,:), ...
                 'FaceAlpha',0.0, 'EdgeColor',ps.colors(2,:), ...
                 'LineStyle',ps.lines(2), 'LineWidth',ps.linewidths(2), ...
                 'DisplayName','$\mathcal X_{data}$');
        end

        xlabel(sprintf('x_{%d}',dim(1)));
        ylabel(sprintf('x_{%d}',dim(2)));
        title(sprintf('k = %d',k-1));
        grid on; box on;
        if exist('axx','var') && numel(axx)>=i_dim && ~isempty(axx{i_dim})
            axis(axx{i_dim});
        end
    end
    lg = legend('Location','best','Interpreter','latex','AutoUpdate','off');
    set(lg.BoxFace,'ColorType','truecoloralpha', ...
        'ColorData',uint8(255*[1;1;1;.9]));
    sgtitle(ps.name);
end

fprintf('[visualizeDDRA]  ✅  Visualization complete.\n');
end
