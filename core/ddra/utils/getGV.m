function Vmatzono = getGV(lookup)
    V = lookup.V;
    totalsamples = lookup.totalsamples; 
    n = lookup.sys.dims.n;

    g = size(V.G,2);  % number of generators
    GV = cell(1, g * totalsamples);  % preallocate
        
    index=1;
    for i=1:g
        vec = V.G(:,i);
        GV{index} = [vec, zeros(n,totalsamples-1)];
        for j=1:totalsamples-1
            GV{j+index}= [GV{index+j-1}(:,2:end) GV{index+j-1}(:,1)];
        end
        index = j+index+1;
    end

    Vmatzono= matZonotope(zeros(n, totalsamples),GV);
end
