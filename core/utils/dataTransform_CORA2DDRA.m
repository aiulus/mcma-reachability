function [x_all, utraj_all] = dataTransform_CORA2DDRA(testSuite)
% DATA-TRANSFORM  Convert a CORA testSuite to the flat matrices expected by
%                 the DDRA routines (full-state output, y = x).
%
%   x_all    ∈ ℝ^{n_y × K·s·T_k}
%   utraj_all∈ ℝ^{n_u × K·s·T_k}

% ---------- basic sizes --------------------------------------------------
K                 = numel(testSuite);                    % # nominal curves
[T_k, n_y, s]     = size(testSuite{1}.y);                % horizon, outputs, MC runs

% ---------- determine input dimension & orientation ---------------------
u0        = testSuite{1}.u;
if ndims(u0) == 2
    u0 = reshape(u0, size(u0,1), size(u0,2), 1);         % promote to 3-D
end
dimsU = size(u0);                                        % [d1 d2 s]

if dimsU(1) == T_k && dimsU(2) ~= T_k
    n_u         = dimsU(2);                              % layout:  T_k × n_u × s
    permuteFlag = false;
elseif dimsU(2) == T_k && dimsU(1) ~= T_k
    n_u         = dimsU(1);                              % layout:  n_u × T_k × s
    permuteFlag = true;                                  % need permute later
else                                                     % both (or neither) equal T_k
    % Pick the smaller one as n_u and assume it is input-first
    [n_u,~]     = min(dimsU(1:2));
    permuteFlag = true;
end
% ------------------------------------------------------------------------

% ---------- pre-allocate output matrices --------------------------------
totTraj  = K * s;
cols     = totTraj * T_k;

x_all     = zeros(n_y, cols);
utraj_all = zeros(n_u, cols);
% ------------------------------------------------------------------------

trajIdx = 0;
for i = 1:K
    y_full = testSuite{i}.y;                             % (T_k × n_y × s)
    u_full = testSuite{i}.u;                             % raw input

    % --- normalise u_full to 3-D, time-first layout ---------------------
    if ndims(u_full) == 2
        u_full = reshape(u_full, size(u_full,1), size(u_full,2), 1);
    end
    if permuteFlag
        u_full = permute(u_full, [2 1 3]);               % → T_k × n_u × s
    end
    % --------------------------------------------------------------------

    for j = 1:s
        trajIdx = trajIdx + 1;
        lb = (trajIdx-1)*T_k + 1;
        ub =  trajIdx   *T_k;

        % outputs
        x_all(:, lb:ub) = squeeze(y_full(:,:,j)).';      % n_y × T_k

        % inputs
        u_slice = squeeze(u_full(:,:,j)).';              % n_u × T_k
        utraj_all(:, lb:ub) = u_slice;
    end
end
end
