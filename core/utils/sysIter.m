function [y_next, y_next_model] = sysIter(y_current, y_model_current, ...
    u, u_model, sys_d, W, V)

    % Sample disturbances
    w_point = randPoint(W);     % Original: w_point = randPoint(W);
    v_point = randPoint(V);     % Original: v_point = randPoint(V);

    % Apply system dynamics
    y_next = sys_d.A * y_current + sys_d.B * u + ...
             w_point + v_point - sys_d.A * v_point;

    % Apply same to model trajectory
    y_next_model = sys_d.A * y_model_current + sys_d.B * u_model + ...
                   w_point + v_point - sys_d.A * v_point;
end
