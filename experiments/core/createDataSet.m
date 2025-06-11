function [data, testSuite] = createDataSet(sys, initpoints, steps, x_all, utraj_all)
    data = aux_DDRA(sys, initpoints, steps, x_all, utraj_all);
    testSuite = aux_CC(data, sys.CORA.dt);
end

function data = aux_DDRA(sys, initpoints, steps, x_all, utraj_all)
%Inputs:
%  sys           -   system struct returned by systemsDDRA()
%  initpoints    -  #(different trajectories simulated in parallel)
%  steps         -  (shared) length of each trajectory
%  x_all         - (initpoints * n) x (steps + 1) returned by getDataDDRA
%  utraj_all     - (initpoints * m) x steps matrix returned by getDataDDRA
% 
% Output:
%  D    =   struct() with fields:
%       .n
%       .m
%       .Ntraj   = initpoins
%       .T       = steps
%       .X0      = (1 x Ntraj) cell array of initial states
%       .Utrajs  = (1 x Ntraj) cell array ofinput sequences
%       .Xtrajs  = (1 x Ntraj) cell array of state sequences
%
% Example:
%  sys = systemsDDRA('chain_of_integrators', 0.05, 5);
% [X0, U, W, Wmatzono] = initialSetupDDRA(sys, 1, 120, ...);
% [x, utraj] = getDataDDRA(sys, 1, 120, X0, U, W);
% D = flatData2TestSuite(sys, 1, 120, x, utraj);

    n = sys.dims.n;
    m = sys.dims.m;
    Ntraj = initpoints;
    T = steps;

    % Preallocate cell arrays
    X0cells   = cell(1, Ntraj);
    Ucells    = cell(1, Ntraj);
    Xcells    = cell(1, Ntraj);

    for i = 1:Ntraj
        % Extract the i-th trajectory from the stacked matrices
        rowStart_x  = (i-1)*n + 1;
        rowEnd_x    = i*n;
        rowStart_u  = (i-1)*m + 1;
        rowEnd_u    = i*m;

        % x_all(rowStart_x:rowEnd_x, :) is n×(T+1) 
        Xcell = x_all(rowStart_x:rowEnd_x, :);   

        % utraj_all(rowStart_u:rowEnd_u, :) is m×T
        Ucell = utraj_all(rowStart_u:rowEnd_u, :);

        % Initial state
        X0cell = Xcell(:, 1);

        % Save into output
        X0cells{i}  = X0cell;        % n×1
        Xcells{i}   = Xcell;         % n×(T+1)
        Ucells{i}   = Ucell;         % m×T
    end

    % Build the commonData struct
    data = struct();
    data.n      = n;
    data.m      = m;
    data.Ntraj  = Ntraj;
    data.T      = T;
    data.X0     = X0cells;  % 1×Ntraj cell, each n×1
    data.Utrajs = Ucells;   % 1×Ntraj cell, each m×T
    data.Xtrajs = Xcells;
end

function testSuite = aux_CC(data, dt)
% Inputs:
%   data must have fields:
%     .n       (state dimension)
%     .m       (input dimension)
%     .Ntraj   (number of separate trajectories)
%     .T       (number of steps per trajectory)
%     .X0      1×Ntraj cell array, each n×1
%     .Utrajs  1×Ntraj cell array, each m×T  (inputs for k=0…T-1)
%     .Xtrajs  1×Ntraj cell array, each n×(T+1)  (states for k=0…T)
%
% Output:
%   testSuite is a 1×Ntraj cell array. Each testSuite{i} has fields:
%       .u  (m×T)   — exactly commonData.Utrajs{i}
%       .y  (n×T×1) — we store x(k) for k=0…T−1 (third dimension = replicates =1)

    Ntraj = data.Ntraj;
    T     = data.T;
    n     = data.n;
    m     = data.m;

    testSuite = cell(1, Ntraj);

    for i = 1:Ntraj
        Yi = data.Xtrajs{i}(:, 1:T);
        Ui = data.Utrajs{i};
        x0 = data.X0{i};
        tc = testCase(Yi, Ui, x0, dt);
        testSuite{i} = tc;
    end
end

