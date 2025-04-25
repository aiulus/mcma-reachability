function y = u2y(u, sys)
    n = size(sys.discrete.A, 1);
    y = inv(eye(n)-sys.discrete.A)*sys.discrete.B*u; 
end

