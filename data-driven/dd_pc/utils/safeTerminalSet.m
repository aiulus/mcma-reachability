%% PLACEHOLDER - will be replaced by safe mRCI set
function S_f = safeTerminalSet(sys, spread_factor, psi_M)
    % Inputs:
    %   sys - system with field `dims.n` (output dim) and `dims.m` (input dim)
    %   spread_factor - scalar controlling bounding set width
    %   psi_M - Number of generators

    n = sys.dims.n;  % output dimension

    % Centered at zero, with generators spread across dimensions
    Y0 = zonotope(zeros(n,1), spread_factor * eye(n, psi_M));

    S_f = Y0;
end

function Z_prod = cartProdZonotope(Z1, Z2)
    % Cartesian product of zonotopes: concatenate centers and generators
    c_combined = [Z1.center; Z2.center];
    G_combined = blkdiag(Z1.generators, Z2.generators);
    Z_prod = zonotope(c_combined, G_combined);
end
