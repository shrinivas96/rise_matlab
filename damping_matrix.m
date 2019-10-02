function D = damping_matrix(x)
    global r b c;
    phi = x(1);
    theta = x(2);
    D = [b/r b/r; 0 c];
end