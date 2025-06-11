function [data, testSuite] = createDataSet(sys, initpoints, steps, x_all, utraj_all)
% createDataSet - Convert flat simulation data into commonData and testCase objects
% Inputs:
%   sys         - system struct from systemsDDRA
%   initpoints  - number of parallel trajectories
%   steps       - trajectory length
%   x_all       - (initpoints*n) x (steps+1) matrix from getDataDDRA
%   utraj_all   - (initpoints*m) x steps matrix from getDataDDRA
% Outputs:
%   data        - struct with fields:
%                   .n, .m, .Ntraj, .T, .X0, .Utrajs, .Xtrajs
%   testSuite   - 1×Ntraj cell array of testCase objects

    % Build commonData
    data = aux_DDRA(sys, initpoints, steps, x_all, utraj_all);
    % Build testSuite of testCase objects
    testSuite = aux_CC(data, sys.CORA.dt);
end

function data = aux_DDRA(sys, initpoints, steps, x_all, utraj_all)
    n = sys.dims.n;
    m = sys.dims.m;
    Ntraj = initpoints;
    T = steps;

    % Preallocate arrays
    X0cells   = cell(1, Ntraj);
    Ucells    = cell(1, Ntraj);
    Xcells    = cell(1, Ntraj);

    for i = 1:Ntraj
        rowStart_x = (i-1)*n + 1;
        rowEnd_x   = i*n;
        rowStart_u = (i-1)*m + 1;
        rowEnd_u   = i*m;
        % Extract data
        Xcell = x_all(rowStart_x:rowEnd_x, :);   % n×(T+1)
        Ucell = utraj_all(rowStart_u:rowEnd_u, :);% m×T
        X0cell = Xcell(:,1);                     % n×1
        X0cells{i} = X0cell;
        Ucells{i}  = Ucell;
        Xcells{i}  = Xcell;
    end

    data = struct();
    data.n      = n;
    data.m      = m;
    data.Ntraj  = Ntraj;
    data.T      = T;
    data.X0     = X0cells;
    data.Utrajs = Ucells;
    data.Xtrajs = Xcells;
end

function testSuite = aux_CC(data, dt)
% aux_CC - build testCase objects with correct dimensions for y and u
    Ntraj = data.Ntraj;
    T     = data.T;
    n     = data.n;
    m     = data.m;
    testSuite = cell(1, Ntraj);
    for i = 1:Ntraj
        % Extract raw trajectories
        Xfull = data.Xtrajs{i};    % n×(T+1)
        Uraw  = data.Utrajs{i};    % m×T
        x0    = data.X0{i};        % n×1
        % Build y: time×outputs×replicates = T×n×1
        Y = permute(Xfull(:,1:T), [2,1]);   % T×n
        Y = reshape(Y, [T, n, 1]);
        % Build u: time×inputs×replicates = T×m×1
        U = permute(Uraw, [2,1]);         % T×m
        U = reshape(U, [T, m, 1]);
        % Create testCase: (y, u, x0, dt)
        tc = testCase(Y, U, x0, dt);
        testSuite{i} = tc;
    end
end
