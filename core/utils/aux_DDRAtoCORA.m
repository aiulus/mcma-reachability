function testSuite = aux_DDRAtoCORA(data, dt)
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

