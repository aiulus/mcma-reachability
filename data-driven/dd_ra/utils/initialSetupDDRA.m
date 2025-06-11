function [X0, U, W, Wmatzono] = initialSetupDDRA(sys, initpoints, steps, ...
    X0_center, X0_spread, U_center, U_spread, W_center, W_spread)

    n = sys.dims.n;
    m = sys.dims.m;
    totalsamples = initpoints*steps;

    X0 = zonotope(X0_center, X0_spread*diag(ones(n, 1)));
    U = zonotope(U_center*ones(m, 1), U_spread*diag(ones(m, 1)));    
    %W = zonotope(W_center*ones(n, 1), W_spread*ones(n, 1));
    W = zonotope(W_center.*ones(n, 1), W_spread.*ones(n, 1));
    % Preallocation for GW
    numGens = size(W.generators,2);
    GW = cell(1, numGens * totalsamples);

    %Construct matrix zonotpe \mathcal{M}_w
    index=1;
    for i=1:numGens
        vec=W.Z(:,i+1);
        GW{index}= [ vec,zeros(n,totalsamples-1)];
        for j=1:totalsamples-1
            GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
        end
        index = j+index+1;
    end
    Wmatzono= matZonotope(zeros(n,totalsamples),GW);
end

