function M_sigma = estimateAB_output(sys, Y_0T, Y_1T, U_0T, M_v, M_w, M_Av)
    %% TODO: Ensure that the pseudo-inverse is well-conditioned 
    Hrank = rank([Y_0T; U_0T]);
    n = size(Y_0T, 1); m = size(U_0T, 1);
    if Hrank < (n + m)
        warning("The matrix [U_- \\ Y_-] doesn't have full row rank: \n " + ...
            "rank([U_- \\ Y_-]) = %d < %d = n + m. \n" + ...
            "The right-inverse for [U_- \\ Y_-] may not be reliable. \n" + ...
            "Consider using a persistently exciting input signal of order n+1 if dealing with noise-free measurements.", Hrank, n + m);
    end
    M_sigma = (Y_1T - M_w - M_v + M_Av)*pinv([Y_0T; U_0T]);
    aux_validate_Msigma(sys, M_sigma);
end

function aux_validate_Msigma(sys, M_sigma)
    intAB = intervalMatrix(M_sigma);
    
    lower_bound = intAB.inf;
    upper_bound = intAB.sup;
    
    true_system = [sys.A, sys.B];
    
    is_contained = all(lower_bound <= true_system, 'all') && all(true_system <= upper_bound, 'all');
    
    if ~is_contained
        warning('Validation failed: M_sigma does not contain the true [A, B].');
    else
        disp('Validation successful: M_sigma contains the true [A, B].');
    end
end

%function aux_validate_Msigma_orig(sys, M_sigma)
%    intAB11 = intervalMatrix(M_sigma);
%    intAB1 = intAB11.int;
%    intAB1.sup >= [sys.A,sys.B]
%    intAB1.inf <= [sys.A,sys.B]
%end