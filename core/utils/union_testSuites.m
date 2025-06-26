function complete_testSuite = union_testSuites(varargin)
    % UNION_TESTSUITES  Concatenate multiple testSuite cell arrays along the replicate axis
    %   complete_testSuite = union_testSuites(ts1, ts2, ..., tsN)
    %   Inputs:
    %     tsi - 1×M cell array of testCase objects, all with identical
    %           .y (T×ny×si), .u (T×nu×si), .initialState (ny×1×si)
    %   Output:
    %     complete_testSuite - 1×M cell array of testCase, where each element's
    %           .y,.u,.initialState,.x are concatenated along the 3rd dim: si = sum(si_i)
    
    % Number of test sets
    numSets = numel(varargin);
    assert(numSets>0, 'Provide at least one testSuite.');
    
    % Determine cell-array size
    M = numel(varargin{1});
    %for k=2:numSets
    %    assert(numel(varargin{k})==M, 'All testSuite inputs must have same length.');
    %end
    
    % Preallocate
    complete_testSuite = cell(1, M);
    
    % Loop over each testCase position
    for i=1:M
        % Start with first
        tc0 = varargin{1}{i};
        Y = tc0.y;       % T×ny×s1
        U = tc0.u;       % T×nu×s1
        X0 = tc0.initialState; % ny×1×s1
        % Also x if present
        if ~isempty(tc0.x)
            Xfull = tc0.x; % ny×(T+1)×s1
        else
            Xfull = [];
        end
        % Concat all others
        for k=2:numSets
            tc = varargin{k}{i};
            Y = cat(3, Y, tc.y);
            U = cat(3, U, tc.u);
            X0 = cat(3, X0, tc.initialState);
            if ~isempty(Xfull)
                Xfull = cat(3, Xfull, tc.x);
            end
        end
        % Rebuild testCase
        tc_new = testCase(Y, U, X0, tc0.sampleTime, tc0.name);
        % Attach x if it was used
        if ~isempty(Xfull)
            tc_new.x = Xfull; 
        end
        complete_testSuite{i} = tc_new;
    end
end
