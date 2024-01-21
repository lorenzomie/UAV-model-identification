function [A,B,C,D] = Model(Xu,Xq,Mu,Mq,X_delta,M_delta,Ts)
% ODE function for computing state-space matrices as functions of parameters
g = 9.81;
A = [Xu Xq -g;
    Mu  Mq 0;
    0 1 0];
B = [X_delta;
    M_delta;
    0];
C = [0 1 0;
     Xu Xq 0]; 

D = [0;X_delta];

end