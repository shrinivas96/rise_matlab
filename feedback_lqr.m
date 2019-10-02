function u = feedback_lqr(x)
    K = [-0.1414    0.4550    0.0904    0.2922];          % for Q = 0.00002*eye(4), R = 0.001;
    u = -K*x;
end