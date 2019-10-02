function A = A_matrix(x)
    
    % global alpha a1 a2 Mlr;
    global alpha Mlr a1 a2;
    phi = x(1);
    theta = x(2);
    % phi_dot = x(3);
    % theta_dot= x(4);
    
    A = [a1 a1+Mlr*cos(theta + alpha); a1+Mlr*cos(theta+alpha) a1+2*Mlr*cos(theta+alpha)+a2];
end