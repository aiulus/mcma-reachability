function noise_catchall = setupInitialNoises(sys, w_spread, v_spread, totalsamples)
    sys_d = sys.discrete;
    n = size(sys_d.A, 1);

    %% PROCESS NOISE
    W = zonotope([zeros(n,1),w_spread*ones(n,1)]);
    
    for i=1:size(W.generators,2)
        vec=W.Z(:,i+1);
        GW{i}= [ vec,zeros(n,totalsamples-1)];
        for j=1:totalsamples-1
            GW{j+i}= [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
        end
    end

    GmatW = cat(3, GW{:});
    WmatZono = matZonotope(zeros(n,totalsamples), GmatW);

    process_noise = struct( ...
        'W', W, ...
        'GW', GW, ...
        'GmatW', GmatW, ...
        'WmatZono', WmatZono ...
        );
    
    %% PROCESS NOISE
    V = zonotope([zeros(n,1),v_spread*ones(n,1)]);
    CV = zeros(n,totalsamples);
    for i=1:size(V.generators,2)
        vec=V.Z(:,i+1);
        GV{i}= [ vec,zeros(n,totalsamples-1)];
        for j=1:totalsamples-1
            GV{j+i}= [GV{i+j-1}(:,2:end) GV{i+j-1}(:,1)];
        end
    end
    GmatV = cat(3, GV{:});
    VmatZono = matZonotope(CV, GmatV);

    msmt_noise = struct( ...
        'V', V, ...
        'GV', GV, ...
        'GmatV', GmatV, ...
        'VmatZono', VmatZono ...
        );
    
    AV = sys_d.A*V;
    VAmatZono = sys_d.A*VmatZono;

    Av = struct( ...
        'AV', AV, ...
        'VAmatZono', VAmatZono ...
        );

    noise_catchall = struct( ...
        'process', process_noise, ...
        'msmt', msmt_noise, ...
        'Av', Av ...
        );
end

