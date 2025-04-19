function u_l = getRandomInput(sys, N, rand_mode, constr_scale)
    m = sys.dims.m;

    lb = sys.intc.inf;
    ub = sys.intc.sup;

    lb = lb*constr_scale; ub = ub*constr_scale;
    
    switch rand_mode
        case 'prbs'
            u_l = idinput([m, N], 'prbs', [0, 1], [-1,1]);
        case 'sinusoid'
            u_l = customSinusoid(lookup.sys, lookup.T_sim, constr_scale, 1e+8);
         case 'sinusoidal_sweep'
            u_l = sinusoidal_sweep(lookup, lb, ub, m, N);       
        case 'uniform'
            u_l = uniform(lb, ub, m, N);
        case 'custom_uniform'
            u_l = customUniform(constr_scale, lb, ub, m, N);
        case 'controlled_random'
            u_l = controlled_random(lookup, lb, ub, m, N);
        case 'white_noise'
            u_l = white_noise(lb, ub, m, N);
    end
end




