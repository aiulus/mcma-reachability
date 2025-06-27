function visualizeReachsets(X0, X_model, X_data, R_cc, projectedDims, numberofplots)
    % visualizeAlinearDT: Visualizes model-based vs data-driven reachable sets.
    %
    % Inputs:
    %   - X0            : initial zonotope
    %   - X_model       : cell array of true reachable sets (R_true)
    %   - X_data        : cell array of data-driven reachable sets (R_ddra)
    %   - R_cc          : cell array of conformance checking sets (R_cc)
    %   - projectedDims : cell array of 2D projection indices (e.g., {[1 2]})
    %   - numberofplots : number of propagation steps to visualize
    
    save_figure = true; % Toggle saving the final plot

    for plotRun = 1:length(projectedDims)
        
        fig_handle = figure('Name', 'Reachable Sets Evolution', 'Renderer', 'painters', 'Position', [10 10 900 700]);
        dims_to_plot = projectedDims{plotRun};

        % --- Auto-calculate axis limits based on all sets ---
        min_x = inf; max_x = -inf;
        min_y = inf; max_y = -inf;
        
        % Combine all sets into a single cell array for easier iteration
        all_sets = [X_model, X_data, R_cc, {X0}];
        
        for i = 1:length(all_sets)
            if ~isempty(all_sets{i}) && iscell(all_sets{i})
                % Handle cases where an element is a cell itself
                current_zono_cell = all_sets{i};
                for j = 1:length(current_zono_cell)
                    if ~isempty(current_zono_cell{j})
                        pZono = project(current_zono_cell{j}, dims_to_plot);
                        intHull = interval(pZono);
                        min_x = min(min_x, infimum(intHull(1)));
                        max_x = max(max_x, supremum(intHull(1)));
                        min_y = min(min_y, infimum(intHull(2)));
                        max_y = max(max_y, supremum(intHull(2)));
                    end
                end
            elseif ~isempty(all_sets{i})
                % Handle single zonotope objects
                pZono = project(all_sets{i}, dims_to_plot);
                intHull = interval(pZono);
                min_x = min(min_x, infimum(intHull(1)));
                max_x = max(max_x, supremum(intHull(1)));
                min_y = min(min_y, infimum(intHull(2)));
                max_y = max(max_y, supremum(intHull(2)));
            end
        end
        
        % Add 10% padding to the axes
        x_range = max_x - min_x;
        y_range = max_y - min_y;
        padding_x = x_range * 0.1;
        padding_y = y_range * 0.1;
        
        axis([min_x - padding_x, max_x + padding_x, min_y - padding_y, max_y + padding_y]);
        hold on; grid on;

        % --- Plotting with color-coding and thicker lines ---
        
        % Plot reachable sets from the true model (R_true)
        for iSet = 1:numberofplots
            handleModel = plot(X_model{iSet}, dims_to_plot, ...
                'Filled', true, 'FaceColor', [0.2 0.2 1], 'FaceAlpha', 0.1, 'EdgeColor', 'b', 'LineWidth', 2);
        end

        % Plot reachable sets from DDRA data (R_ddra)
        for iSet = 1:numberofplots
            handleData = plot(X_data{iSet}, dims_to_plot, 'r-', 'LineWidth', 2);
        end

        % Plot reachable sets from Conformance Checking (R_cc)
        for iSet = 1:length(R_cc)
            handleCC = plot(R_cc{iSet}, dims_to_plot, 'g-', 'LineWidth', 2);
        end
        
        % Plot initial set on top
        handleX0 = plot(X0, dims_to_plot, 'k-', 'LineWidth', 2.5);

        % --- Final Touches: Labels and Legend ---
        xlabel(['x_{', num2str(dims_to_plot(1)), '}']);
        ylabel(['x_{', num2str(dims_to_plot(2)), '}']);
        title('Evolution of True, DDRA, and CC Reachable Sets');
        
        % Prevent legend warnings and create a comprehensive legend
        warOrig = warning; warning('off','all');
        legend([handleX0, handleModel, handleData, handleCC], ...
            'Initial Set', 'True Reach-Set (R_{true})', 'DDRA Reach-Set (R_{ddra})', 'CC Reach-Set (R_{cc})', ...
            'Location', 'northwest', 'FontSize', 12);
        warning(warOrig);
        
        % Style axes
        ax = gca;
        ax.FontSize = 14;
        box on;
    end
    
    fprintf('[visualizeAlinearDT] ✅ Visualization complete.\n');
    
    if save_figure
        outputDir = '../outputs/figures';
        if ~exist(outputDir, 'dir')
           mkdir(outputDir)
        end
        saveas(fig_handle, fullfile(outputDir, 'reachsets_evolution.png'));
        fprintf('Saved figure to %s\n', fullfile(outputDir, 'reachsets_evolution.png'));
    end
end
