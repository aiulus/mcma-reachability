function X_data = ddra_reachsets_poly(AB, params, options)    
    intAB11 = intervalMatrix(AB);
    intAB1 = intAB11.int;

    N = options.tFinal;
    
    X_data = cell(N+1,1);
    X_data{1} = params.R0;
    for i = 1:N
    
    % monomial of intervals of reachable sets and input
    X_z1 = interval(X_data{i});
    U_int = interval(params.U);
    
    cardint =zonotope([interval([1]);X_z1;X_z1.*X_z1;X_z1(1)*X_z1(2);U_int;U_int.*U_int;U_int(1)*U_int(2);X_z1.*U_int;X_z1(1)*U_int(2);X_z1(2)*U_int(1)]);
    X_data{i+1} =AB *cardint + options.W;
    X_data{i+1,1}=reduce(X_data{i+1,1},'girard',options.zonotopeOrder);
    end
    % Reachability Analysis ---------------------------------------------------
    
end