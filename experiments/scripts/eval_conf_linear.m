system_type = 'chain_of_integrators';
dim_list = 2:2:10;
noise_scales = [0.0, 0.1, 0.5, 1.0, 1.5, 2];
plot_toggle = false;

results = struct();

for d = 1:length(dim_list)
    for ns = 1:length(noise_scales)
        dim = dim_list(d);
        noise = noise_scales(ns);

        fprintf('Evaluating CORA Conformance for dim = %d, noise = %.2f\n', dim, noise);
        out = func_conformance_linear(system_type, dim, noise, plot_toggle);
        
        results.time(d, ns) = out.time;
        % TODO: Compra reachset sizes
    end
end
