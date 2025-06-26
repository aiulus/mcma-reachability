function AB = estimate_AB_poly(X_0T, X_1T, U_0T, WmatZ)
    %compute monomials of trajectories
    totalsamples = size(X_0T, 2);
    X_2 = X_0T .*X_0T;
    X1X2 = X_0T(1,:) .* X_0T(2,:);
    U_2 = U_0T .*U_0T;
    U1U2 = U_0T(1,:) .* U_0T(2,:);
    XU = X_0T .*U_0T;
    X1U2X2U1 = X_0T .* U_0T([2 1],:);
    data=[ones(1,totalsamples);X_0T;X_2;X1X2;U_0T;U_2;U1U2;XU;X1U2X2U1];
    
    % set of A B that is consistent with the data
    %AB = (X_1T - WmatZ)*pinv(data);
    % Temporarily setting noise to zero
    AB  = X_1T * pinv(data);
end
% ------------------------------ END OF CODE ------------------------------