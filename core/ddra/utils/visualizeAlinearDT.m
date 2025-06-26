function visualizeAlinearDT(X0, X_model, X_data, projectedDims, axx, numberofplots)
    % visualizeAlinearDT: Visualizes model-based vs data-driven reachable sets.
    %
    % Inputs:
    %   - X0            : initial zonotope
    %   - X_model       : cell array of reachable sets from model
    %   - X_data        : cell array of reachable sets from data
    %   - projectedDims : cell array of 2D projection indices (e.g., {[1 2], [3 4]})
    %   - axx           : axis bounds for each projection (cell array)
    %   - numberofplots : number of propagation steps to visualize
    %

    for plotRun = 1:length(projectedDims)
        
        figure('Renderer', 'painters', 'Position', [10 10 700 900]);

        % Plot initial set
        handleX0 = plot(X0, projectedDims{plotRun}, 'k-', 'LineWidth', 2);
        hold on;

        % Plot reachable sets from model
        for iSet = 2:numberofplots
            handleModel = plot(X_model{iSet}, projectedDims{plotRun}, ...
                'b', 'Filled', true, 'FaceColor', [.8 .8 .8], 'EdgeColor', 'b');
        end

        % Plot reachable sets from data
        for iSet = 2:numberofplots
            handleData = plot(X_data{iSet}, projectedDims{plotRun}, 'r');
        end

        % Label axes
        xlabel(['x_{', num2str(projectedDims{plotRun}(1)), '}']);
        ylabel(['x_{', num2str(projectedDims{plotRun}(2)), '}']);

        % Optional axis setting (if defined)
        if nargin >= 5 && length(axx) >= plotRun
            axis(axx{plotRun});
        end

        % Handle legends without warning spam
        warOrig = warning; warning('off','all');
        legend([handleX0, handleModel, handleData], ...
            'Initial Set', 'Set from Model', 'Set from Data', 'Location', 'northwest');
        warning(warOrig);

        % Style axes
        ax = gca;
        ax.FontSize = 22;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3) - 0.01;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
    end

    fprintf('[visualizeAlinearDT] ✅ Visualization complete.\n');
end
