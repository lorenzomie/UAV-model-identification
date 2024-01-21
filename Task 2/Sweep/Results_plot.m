clc; clear; close all
load('MONTECARLO_SW')
set(0, 'DefaultLineLineWidth', 1.5); set(0,'defaultAxesFontSize',16);

%% Plot traces of covariance matrices of Sine Sweep

error_mean = mean(MONTECARLO_SW.trace);
error_variance = var(MONTECARLO_SW.trace);

semilogy(MONTECARLO_SW.trace)
title('Covariance matrix Trace')
xlabel('Samples')
grid on

%% Bode plot firs Montecarlo value (bad estimation)

param_est = cell2mat(MONTECARLO_SW.param_estimated(1));

Xu = param_est(1);
Xq = param_est(2);
Mu = param_est(3);
Mq = param_est(4);
Xd = param_est(5);
Md = param_est(6);

g = 9.81;
A = [Xu Xq -g;
    Mu  Mq 0;
    0 1 0];
B = [Xd;
    Md;
    0];
C = [0 1 0;
     Xu Xq 0]; 

D = [0;Xd];

figure
usys = ss(A,B,C,D);
bode(usys,'b');
legend(Location="best")
grid
set(findall(gcf,'type','line'),'linewidth',1.5)
