function [X_minus, X_plus, U_plus] = shift_trajs(T_k, x, utraj)
    n = size(x, 1); m = size(utraj, 1);
    totalsamples = size(x, 2);
    totalsamples = (totalsamples / T_k)* (T_k - 1);
    
    X_minus = zeros(n, totalsamples);     
    X_plus = zeros(n, totalsamples);
    U_plus = zeros(m, totalsamples);

    for j=1:T_k
        lower = 1 + (j - 1)*(T_k - 1);
        upper = j * (T_k - 1);
        minus_lower = (j - 1) * T_k + 1; minus_upper = (j - 1) * T_k + (T_k - 1);
        plus_lower = (j - 1) * T_k + 2; plus_upper = j * T_k;
        X_minus(:, lower:upper) = x(:, minus_lower:minus_upper);
        X_plus(:, lower:upper) = x(:, plus_lower:plus_upper);
        U_plus(:, lower:upper) = utraj(:, plus_lower:plus_upper);
    end
end

