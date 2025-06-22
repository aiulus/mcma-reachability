function [x,y] = aux_testSuite2regress(testSuite, p)
    % Corrected aux_testSuite2regress function
    % x = [y_1(k-p) ... y_n(k-p) ... y_1(k-1) ... y_n(k-1) ...
    %      u_1(k-p) ... u_n(k-p) ... u_1(k) ... u_n(k)]
    
    if isempty(testSuite)
        x = [];
        y = [];
        return;
    end
    
    % Extract dimensions from the first test case
    y_first = testSuite{1}.y;
    u_first = testSuite{1}.u;
    dim_y = size(y_first, 2);
    dim_u = size(u_first, 2);
    n_k_first = size(y_first,1);
    n_s_first = size(y_first,3);
    
    % Calculate expected number of features and total size
    n_features = dim_y*p + dim_u*(p+1);
    total_size = length(testSuite) * (n_k_first - p) * n_s_first;
    
    x = zeros(total_size, n_features);
    y = zeros(total_size, dim_y);
    index = 1;
    
    for m = 1:length(testSuite)
        y_m = testSuite{m}.y;
        u_m = testSuite{m}.u;
        n_s = size(y_m, 3); % number of samples in this test case
    
        for k = p+1:size(y_m,1)
            
            % Correctly reshape y_features per sample
            y_features = reshape(permute(y_m(k-p:k-1,:,:), [2 1 3]), 1, [], n_s);
            
            % Correctly reshape u_features per sample
            % This handles cases where u_m is 2D (same u for all samples) or 3D
            if ndims(u_m) == 2
                u_features_single = reshape(permute(u_m(k-p:k,:), [2 1]), 1, []);
                u_features = repmat(u_features_single, [1, 1, n_s]);
            else
                u_features = reshape(permute(u_m(k-p:k,:,:), [2 1 3]), 1, [], n_s);
            end
            
            % Concatenate features
            x_k = [y_features, u_features];
    
            y_k = y_m(k,:,:);
    
            % Assign to the main matrices
            end_index = index + n_s - 1;
            x(index:end_index, :) = squeeze(x_k)';
            y(index:end_index, :) = squeeze(y_k)';
            index = end_index + 1;
        end
    end
end