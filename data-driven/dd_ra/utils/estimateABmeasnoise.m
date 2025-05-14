function [AB_av, AB_cmz] = estimateABmeasnoise(U_full, X_0T, X_1T, lookup)
    %% Parameter extraction
    n = lookup.sys.dims.m;
    Wmatzono = lookup.Wmatzono; 
    Vmatzono = lookup.Vmatzono;
    AVmatzono = lookup.AVmatzono;
    totalsamples = lookup.totalsamples;
    
    % Minkowski sum
    mink= -1*Wmatzono +  -1*Vmatzono + AVmatzono;

    % Compute A_tildeNsigma and B_tildeNsigma
    basis=null([X_0T;U_full]);

    % Pre-allocate Acon
    g = numel(mink.G);
    Acon = cell(1, g);

    for i=1:length(mink.generator)
        Acon{i}=(mink.generator{i})*basis;
    end

    Bcon = (X_1T - zeros(n, totalsamples))*basis  ;
    
    % AB with AV assumption (tildeMsigma)
    AB_av = (X_1T + -1*Wmatzono + -1* Vmatzono + AVmatzono)*pinv([X_0T; U_full]);    
    
    % constrained matrix zonotope AB with AV assumption (tildeNsigma)
    AB_cmz = conMatZonotope(AB_av.c, AB_av.G, Acon,Bcon);
end

