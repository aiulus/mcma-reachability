function [sys, R0, U, p_true] = getSystemConf(sys_type, dt, dim)
% SYSTEMSZONODDSF 
% - takes in a string, returns both continuous and discrete-time system matrices 
% - main responsibility: defining various well-known, well-studied systems
%
%   Inputs:
%       SYS_TYPE - A string specifying the type of system to initialize.
%                  Available options include:
%                    - 'inverted_pendulum'
%                    - 'damper'
%                    - 'quadrotor'
%                    - 'dc_motor'
%                    - 'cruise_control'
%                    - 'acc'
%                    - 'ballNbeam'
%                    - 'double_pendulum'
%
%       DT  - Sampling time
%
%   Outputs:
%       SYS - A structure containing:
%               - cont: Continuous-time state-space model (A, B, C, D)
%               - discrete: Discrete-time state-space model (Ad, Bd, Cd, Dd)
%
%   See also: c2d

    if nargin < 3 && strcmpi(sys_type,'chain_of_integrators')
        error('Target system dimension must be specified for systype <chain of integrators>!');
    end

    % DEFAULT-- Change individually within switch-case
    spread_factors = struct( ...
        'R0', 1, ...
        'W', 1, ...
        'V', 1, ...
        'U', 1 ...
        );

    switch lower(sys_type)        
        case 'chain_of_integrators'
            theta = []; % No parameters
            n = dim;
            A = diag(ones(n-1, 1), 1);
            B = zeros(n, 1); B(end) = 1;
            C = eye(n);
            D = zeros(n, 1);

        case 'example0'
            theta = []; % No parameters
            A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
            B = ones(5,1);
            C = [1,0,0,0,0];
            D = 0;
        case 'quadrotor'
            theta = []; % Placeholder
            mass = 0.2;
            g = 9.81;
            I = repmat(1e-3, 3, 1);
            n = 12;
            m_input = 4;
            
            A_i = [1, 2, 3, 10, 11, 12, 8, 7];
            A_j = [4, 5, 6, 7, 8, 9, 1, 2];
            A_val = [ones(6, 1); g; -g];
            A = sparse(A_i, A_j, A_val, n, n);
    
            B_i = [9, 4, 5, 6];
            B_j = [1, 2, 3, 4];
            B_val = [1/mass, 1/I(1), 1/I(2), 1/I(3)];
            B = sparse(B_i, B_j, B_val, n, m_input);
    
            indices = [1, 2, 3, 10, 11, 12];
            C = sparse(1:length(indices), indices, 1, length(indices), n);
            D = zeros(6, 4);
    
        case 'damper'
            m = 100;
            k = 1;
            b = 0.2;

            theta = [m, k, b];
    
            A = [0 1; -k/m -b/m];
            B = [0; 1/m];
            C = [1 0];
            D = 0;
    
        case 'inverted_pendulum'
            M = 50;
            m = 2;
            I = 0.6;
            l = 3;
            b = 0.1;
            g = 9.81;

            theta = [M, m, I, l, b];
    
            p = I*(M+m)+M*m*l^2;
    
            A = [0      1              0           0;
                 0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
                 0      0              0           1;
                 0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
    
            B = [     0;
                 (I+m*l^2)/p;
                      0;
                   m*l/p];
    
            C = [1 0 0 0;
                 0 0 1 0];
            D = [0; 0];
    
        case 'dc_motor'
            J = 0.01;
            b = 0.1;
            K = 0.01;
            R = 1;
            L = 0.5;

            theta = [J, b, K, R, L];
    
            A = [-b/J K/J; -K/L -R/L];
            B = [0; 1/L];
            C = [1 0];
            D = 0;
    
        case 'cruise_control'
            mass = 1000;
            damping = 50;

            theta = [mass, damping];
    
            A = -damping / mass;
            B = 1 / mass;
            C = 1;
            D = 0;
    
        case 'acc'
            mc = 1650;

            theta = mc;
    
            A = [0 1; 0 0];
            B = [0; 1/mc];
            C = [1 0];
            D = 0;
    
        case 'ballnbeam'
            m = 0.11;
            R = 0.015;
            d = 0.03;
            L = 1;
            J = 9.99e-6;
            g = 9.8;

            theta = [m, R, d, L, J];
    
            b21 = - (m * g * d) / (L * (m + (J / R^2)));
    
            A = [0 1; 0 0];
            B = [0; b21];
            C = [1 0];
            D = 0;
    
        case 'double_pendulum'
            M = 100;
            M1 = 10;
            M2 = 10;
            L1 = 2;
            L2 = 1;
            g = 9.81;

            theta = [M ,M1, M2, L1, L2];
    
            A = sparse([
                0, 0, 0, 1, 0, 0;
                0, 0, 0, 0, 1, 0;
                0, 0, 0, 0, 0, 1;
                0, -g*(L1*(M1+M2)+L2*M2)/(L1*M), -g*(L2*M2)/(L1*M), 0, 0, 0;
                0, g*(L1*M1*(M + M1 +M2)+L2*M2*(M+M1))/(M*M1*L1^2), ...
                   g*M2*(-L1*M+L2*(M+M1))/(M*M1*L1^2), 0, 0, 0;
                0, -g*M2/(L1*M1), ...
                   g*(L1*(M1+M2)-L2*M2)/(L1*L2*M1), 0, 0, 0
            ]);
    
            B = sparse([
                0;
                0;
                0;
                1/M1;
                -1/(L1*M);
                0
            ]);
    
            C = sparse([
                1, 0, 0, 0, 0, 0;
                1, L1, 0, 0, 0, 0;
                1, L1+L2, L2, 0, 0, 0
            ]);
    
            D = sparse(zeros(3, 1));
    
        otherwise
            error('Unknown system type: %s', sys_type);
    end
    
    % Discretize using c2d
    sys_cont = ss(A, B, C, D);
    sys_disc = c2d(sys_cont, dt);
    dims = struct( ...
        'n', size(A, 1), ...
        'm', size(B, 2), ...
        'p', size(C, 1) ... % MUST EQUAL n IN OUR CURRENT SETTING!
        );

    % Populate the output structure    
    sys = struct( ...
        'cont', struct('A', A, 'B', B, 'C', C, 'D', D), ...
        'discrete', struct('A', sys_disc.A, 'B', sys_disc.B, 'C', sys_disc.C, 'D', sys_disc.D), ...
        'dims', dims, ...
        'systype', sys_type ...
    );

    p_true = theta;
    [R0, U] = build_std_uncertainty_sets(dims, spread_factors);
end

function [R0, U] = build_std_uncertainty_sets(dims, spread_factors)
    dim_x = dims.n; dim_y = dims.p; dim_u = dims.m;

    R0 = zonotope(zeros(dim_x, 1), spread_factors.R0*eye(dim_x));
    V = zonotope(zeros(dim_y, 1), spread_factors.V*eye(dim_y));
    W = zonotope(zeros(dim_x, 1), spread_factors.W*eye(dim_x));
    % TODO: Probably unnecessary to define V, W first if the center vectors
    %       & generators are hard-coded
    U = cartProd(W, V);
end



