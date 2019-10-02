function B = B_matrix(x)
    
    global alpha n1 Mlr;
    % phi = x(1);
    theta = x(2);
    % phi_dot = x(3);
    theta_dot= x(4);
    n1 = -1*(Mlr*theta_dot*sin(theta + alpha));
    B = [0, n1; 0, n1];
end
