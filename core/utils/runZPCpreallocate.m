function [y_t, y_t_model, uPred, uPred_model, ...
          YPred, YPred_model, execTimeZPC, execTimeRMPC, ...
          Rplotall, yt2ref, yt2ref_model, ...
          R_N, R_N_int, y_N, u_N] = ...
          runZPCpreallocate(p, N, maxsteps)

    %% --- Preallocation steps from original function ---
    % From: p = size(y0,1);
    y_t         = zeros(p, maxsteps + 1);      % y_t = zeros(p, maxsteps + 1);
    y_t_model   = zeros(p, maxsteps + 1);      % y_t_model = zeros(p, maxsteps + 1);

    uPred       = zeros(1, maxsteps);          % uPred = zeros(1, maxsteps);
    uPred_model = zeros(1, maxsteps);          % uPred_model = zeros(1, maxsteps);

    YPred       = zeros(p, maxsteps + 1);      % YPred = zeros(p, maxsteps + 1);
    YPred_model = zeros(p, maxsteps + 1);      % YPred_model = zeros(p, maxsteps + 1);

    execTimeZPC  = zeros(1, maxsteps);         % execTimeZPC = zeros(1, maxsteps);
    execTimeRMPC = zeros(1, maxsteps);         % execTimeRMPC = zeros(1, maxsteps);

    Rplotall = cell(1, maxsteps);              % Rplotall = cell(1, maxsteps);

    yt2ref       = zeros(1, maxsteps);         % yt2ref = zeros(1, maxsteps);
    yt2ref_model = zeros(1, maxsteps);         % yt2ref_model = zeros(1, maxsteps);

    % Optional plotting data
    R_N     = cell(1, N+1);                    % R_N = cell(1, N+1);
    R_N_int = cell(1, N+1);                    % R_N_int = cell(1, N+1);
    y_N     = cell(1, N+1);                    % y_N = cell(1, N+1);
    u_N     = cell(1, N);                      % u_N = cell(1, N);
end
