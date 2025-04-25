function [y_t, y_t_model, YPred] = runZPCinitializeY0(y0, y_t, y_t_model, YPred, k)
    % -- From original code: if k == 1 ... --
    y_t(:, k)       = y0;       % Original: y_t(:,k) = y0;
    y_t_model(:, k) = y0;       % Original: y_t_model(:,k) = y0;
    YPred(:, 1)     = y0;       % Original: YPred(:,1) = y0;
end
