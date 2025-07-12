function M_sigma = estimateAB_output(Y_0T, Y_1T, U_full, M_v, M_w, M_Av)
    M_sigma = (Y_1T - M_v - M_w + M_Av)+pinv([Y_0T; U_full]);
end

