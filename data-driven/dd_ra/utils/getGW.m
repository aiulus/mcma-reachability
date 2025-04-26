function Wmatzono = getGW(lookup)
    W = lookup.W;
    totalsamples = lookup.totalsamples; 
    n = lookup.n;

    g = size(W.generators,2);  % number of generators
    GW = cell(1, g * totalsamples);  % preallocate the cell array
    
    index = 1;
    for i = 1:g
        vec = W.Z(:,i+1);
        for j = 0:totalsamples-1
            GW{index} = [zeros(n,j), vec, zeros(n,totalsamples-j-1)];
            index = index + 1;
        end
    end
    
    Wmatzono = matZonotope(zeros(n,totalsamples), GW);

end