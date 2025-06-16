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
        
        %cardint =zonotope([interval([1]);X_z1;X_z1.*X_z1;X_z1(1)*X_z1(2);U_int;U_int.*U_int;U_int(1)*U_int(2);X_z1.*U_int;X_z1(1)*U_int(2);X_z1(2)*U_int(1)]);
        % Inputs
        x1 = X_z1(1);
        x2 = X_z1(2);
        u1 = U_int(1);
        u2 = U_int(2);
        
        % Polynomial features (1st and 2nd order terms)
        features = [ ...
            interval([1]);        % Constant term
            X_z1;                 % x1, x2
            X_z1.^2;              % x1^2, x2^2
            x1 * x2;              % Cross term x1*x2
            U_int;                % u1, u2
            U_int.^2;             % u1^2, u2^2
            u1 * u2;              % Cross term u1*u2
            X_z1 .* U_int;        % x1*u1, x2*u2
            x1 * u2;              % Cross term x1*u2
            x2 * u1               % Cross term x2*u1
        ];
        
        % Construct the zonotope   
        %cardint = zonotope(cardint);
        % Initialize center and generator matrix
        dim = length(features);
        center = zeros(dim, 1);
        G = zeros(dim, dim);  % One generator per dimension
        
        for k = 1:dim
            intv = features(k);  % extract interval
            lb = min(intv);      % lower bound
            ub = max(intv);      % upper bound
            center(k) = (lb + ub) / 2;
            G(k, k) = (ub - lb) / 2;  % single generator per dimension
        end
        
        % Construct the zonotope
        cardint = zonotope(center, G);
    
        X_data{i+1} = AB *cardint + options.W;
        X_data{i+1,1} = reduce(X_data{i+1,1},'girard',options.zonotopeOrder);
    end
    % Reachability Analysis ---------------------------------------------------
    
end