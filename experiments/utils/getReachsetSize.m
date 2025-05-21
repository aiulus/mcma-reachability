function Z_size = getReachsetSize(Z, method)
% computeReachsetSize - Compute a size measure of a zonotope
%
% Inputs:
%   Z      - Zonotope object (assumed to have fields .center, .generators)
%   method - String specifying the type of measure:
%            'volume'       : approximate volume using determinant of GG^T
%            'bbox_volume'  : volume of axis-aligned bounding box
%            'fro_norm'     : Frobenius norm of the generator matrix
%
% Output:
%   size       - scalar size of the reachable set

    G = Z.generators;

    switch lower(method)
        case 'volume'
        % Approximate volume as sqrt(det(G*G^T)) if G*G^T is full-rank
        GG = G*G';
        if rank(GG) == size(GG, 1)
            Z_size = sqrt(det(GG));
        else
            warning(['Squared generator matrix not full-rank! ' ...
                'Volume potentially underestimated.']);
            Z_size = sqrt(det(GG + 1e-6 * eye(Z_size(GG))));
        end

        case 'bbox_volume'
            % Max absolute sum of generators
            max_extent = sum(abs(G), 2);
            Z_size = prod(2 * max_extent);

        case 'fro_norm'
            Z_size = norm(G, 'fro');

        otherwise
            error('Unknown method for reachset size: %s', method);
    end

end

