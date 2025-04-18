function MZ = zonotopeToMatZonotope(Z, dim, T)
% zonotopeToMatZonotope - Creates a matrix zonotope from a scalar zonotope.
%
% Inputs:
%   Z   - Zonotope with center in R^dim and generator matrix in R^{dim x g}
%   dim - Dimension of the state/output
%   T   - Number of time steps (columns of the matrix)
%
% Output:
%   MZ - matZonotope object of dimension (dim x T)

    numGen = size(Z.generators, 2);
    base_center = zeros(dim, T);
    Gcell = cell(1, numGen * T);

    for i = 1:numGen
        vec = Z.Z(:, i+1);
        G = [vec, zeros(dim, T - 1)];

        for j = 1:T
            idx = (i - 1) * T + j;
            Gcell{idx} = circshift(G, [0, j-1]);
        end
    end

    Gmat = cat(3, Gcell{:});
    MZ = matZonotope(base_center, Gmat);
end
