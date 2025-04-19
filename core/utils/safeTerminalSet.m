%% PLACEHOLDER - will be replaced by safe mRCI set
function S_f = safeTerminalSet(sys, U, Y, spread_factor)
    % Inputs:
    %   sys - system with field `dims.n` (output dim) and `dims.m` (input dim)
    %   U   - zonotope representing admissible inputs
    %   Y   - zonotope representing admissible outputs
    %   spread_factor - scalar controlling bounding set width

    %% 1) Extract system dimensions
    n = sys.dims.n;  % output dimension
    m = sys.dims.m;  % input dimension

    %% 2) Define bounding zonotopes U0 and Y0 around origin
    psi_M = 10;  % number of generators

    % Centered at zero, with generators spread across dimensions
    U0 = zonotope(zeros(m,1), spread_factor * eye(m, psi_M));
    Y0 = zonotope(zeros(n,1), spread_factor * eye(n, psi_M));

    %% 3) Compute cross products (U x Y) and (U0 x Y0)
    % Cross product is cartesian product of zonotopes, i.e., concatenation
    UxY  = cartProdZonotope(U, Y);   % Combined safe set
    U0xY0 = cartProdZonotope(U0, Y0); % Bounding invariant template

    %% 4) Compute intersection to form terminal safe set
    S_f = intersect(UxY, U0xY0);  % Compute intersection of zonotopes
end

function Z_prod = cartProdZonotope(Z1, Z2)
    % Cartesian product of zonotopes: concatenate centers and generators
    c_combined = [Z1.center; Z2.center];
    G_combined = blkdiag(Z1.generators, Z2.generators);
    Z_prod = zonotope(c_combined, G_combined);
end
