function sys = setupBoundaryConditions(sys)
% ASSIGNBOUNDARYCONDITIONS Assigns boundary condition sets to the system struct
%   Inputs:
%       sys - system struct from systemsZonoDDSF()
%   Outputs:
%       sys - updated with sys.bcs struct

    sys_type = lower(sys.systype);
    dims = sys.dims;
    n = dims.n;
    m = dims.m;

    switch sys_type
        case 'example0'
            % Output constraints
            y_lb = [-10;2;-10;-10;-10];  
            y_ub = [10;10;10;10;10];     

            % Initial condition
            y0 = [-2;4;3;-2.5;5.5];      % Placeholder
            u0 = 8;

            % Initial state zonotope
            initial_state_spread = 25;    % Placeholder           

            initial_input_spread = 19;           % Placeholder (20 - 1)
            
            
        case 'quadrotor'
            % TODO

        case 'dc_motor'
            % TODO   

        case 'damper'
            % TODO

        case 'inverted_pendulum'
            % TODO    

        case 'cruise_control'
            % TODO

        case 'acc'
            % TODO

        case 'ballnbeam'
            % TODO

        case 'double_pendulum'
            % TODO
            
        otherwise
            error('assignBoundaryConditions: unknown system type %s', sys_type);
    end

    [X0, U, Y, intc] = applySpreadFactors(initial_state_spread, initial_input_spread, y_lb, y_ub, y0, u0, n);

    % Package into struct
    bcs = struct( ...
        'intc', intc, ...
        'y0', y0, ...
        'X0', X0, ...
        'U', U, ...
        'Y', Y ...
    );
    sys.bcs = bcs; % Assign bcs to system
    validate_bcs(sys.bcs, sys.dims); % Validate
end

function [X0, U, Y, intc] = applySpreadFactors(initial_state_spread, initial_input_spread, y_lb, y_ub, y0, u_0, n)
    intc = interval(y_lb, y_ub); % CORA interval object
    X0 = zonotope([y0, initial_state_spread * diag(ones(n,1))]);
    U = zonotope([u_0 - 1, initial_input_spread]);            
    Y = zonotope(intc);
end

function validate_bcs(bcs, dims)
% VALIDATE_BCS Checks the validity of boundary condition definitions

    dim_x = dims.n;
    dim_u = dims.m;

    % Check y0
    assert(isvector(bcs.y0) && length(bcs.y0) == dim_x, ...
        'y0 must be a vector of size %d', dim_x);

    % Check X0 (zonotope center)
    assert(isa(bcs.X0, 'zonotope'), 'X0 must be a CORA zonotope');
    assert(length(bcs.X0.center) == dim_x, ...
        'X0 center must be of length %d', dim_x);
    assert(size(bcs.X0.generators,1) == dim_x, ...
        'X0 generators must have %d rows', dim_x);

    % Check U (input zonotope)
    assert(isa(bcs.U, 'zonotope'), 'U must be a CORA zonotope');
    assert(length(bcs.U.center) == dim_u, ...
        'U center must be of length %d', dim_u);
    assert(size(bcs.U.generators,1) == dim_u, ...
        'U generators must have %d rows', dim_u);

    % Check intc (interval object)
    assert(isa(bcs.intc, 'interval'), 'intc must be a CORA interval');
    assert(length(bcs.intc.inf) == dim_x && length(bcs.intc.sup) == dim_x, ...
        'intc must match system state/output dimension %d', dim_x);
end
