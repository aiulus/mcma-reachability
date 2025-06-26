function [X_0T, X_1T, U_0T, U_1T, U_full] = shiftTrajs_CORA2DDRA(testSuite)
% shiftTrajs_CORA2DDRA
% ---------------------------------------------------------------
% Converts a CORA testSuite into the *shifted* data matrices used
% by DDRA-style identification:
%   X_0T  –  x_k      (dim_y × K·s·(T_k-1))
%   X_1T  –  x_{k+1}
%   U_0T  –  u_k      (dim_u × …)
%   U_1T  –  u_{k+1}
%   U_full–  full input traces (dim_u × K·s·T_k)
%
% Assumption : full observability  (y = x)

    K                 = numel(testSuite);               % #nominal curves
    [T_k, dim_y, s]   = size(testSuite{1}.y);           % horizon, outputs, samples
    T_km1             = T_k - 1;

    % ---- robust input-dimension detection (same logic as dataTransform) ---
    size_u0 = size(testSuite{1}.u);
    if numel(size_u0) == 3            % linearSysDT path   (n_u × T_k × s)
        dim_u = size_u0(1);
    elseif size_u0(1) == T_k          % nonlinear/NARX     (T_k × n_u)
        dim_u = size_u0(2);
    else
        dim_u = size_u0(1);           % rare fallback
    end
    % -----------------------------------------------------------------------

    nTraj         = K * s;                                % total trajectories
    cols_full     = nTraj * T_k;                          % per-step matrices
    cols_shift    = nTraj * T_km1;                        % per-step minus 1

    % -------- pre-allocate --------------------------------------------------
    X_0T  = zeros(dim_y, cols_shift);
    X_1T  = zeros(dim_y, cols_shift);
    U_0T  = zeros(dim_u, cols_shift);
    U_1T  = zeros(dim_u, cols_shift);
    U_full= zeros(dim_u, cols_full);
    % -----------------------------------------------------------------------

    trajIdx = 0;                                          % running counter

    for i = 1:K
        tc       = testSuite{i};
        y_full   = tc.y;                                  % (T_k × dim_y × s)
        u_full   = tc.u;                                  % see below

        % ----- make u_full 3-D with shape (T_k × dim_u × s) ----------------
        if ndims(u_full) == 2
            u_full = repmat(u_full, [1 1 s]);             % duplicate nominal input
        end
        if size(u_full,1) ~= T_k                          % linear case → permute
            u_full = permute(u_full, [2 1 3]);
        end
        % -------------------------------------------------------------------

        for j = 1:s                                        % local samples only
            trajIdx      = trajIdx + 1;

            % global column indices
            lb_full      = (trajIdx-1)*T_k   + 1;
            ub_full      =  trajIdx   *T_k;
            lb_shift     = (trajIdx-1)*T_km1 + 1;
            ub_shift     =  trajIdx   *T_km1;

            % -------- outputs ---------------------------------------------
            y = squeeze(y_full(:,:,j));                    % (T_k × dim_y)
            X_series = y.';                                % dim_y × T_k
            X0 = X_series(:,1:end-1);
            X1 = X_series(:,2:end);
            X_0T(:, lb_shift:ub_shift) = X0;
            X_1T(:, lb_shift:ub_shift) = X1;

            % -------- inputs ----------------------------------------------
            u = squeeze(u_full(:,:,j)).';                  % dim_u × T_k
            U_full(:, lb_full:ub_full) = u;
            U_0T(:, lb_shift:ub_shift) = u(:,1:end-1);
            U_1T(:, lb_shift:ub_shift) = u(:,2:end);
        end
    end
end
