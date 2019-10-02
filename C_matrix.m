function C = C_matrix(x)
    
    global M m g r alpha Mgl;
    phi = x(1);
    theta = x(2);
    c11 = (M+m)*g*r*sin(alpha);
    C = [c11; c11 + Mgl*sin(theta)];
end
