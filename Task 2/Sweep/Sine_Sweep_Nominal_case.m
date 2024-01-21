%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 13/12/2017                                                        %
% Adapted to ANT-X 2DoF by:  Salvatore Meraglia (salvatore.meraglia@polimi.it)%
% Date: 22/12/2022                                                        %
%
% Further modified to include structure three-state identified longitudinal model
% 06/01/23 ML
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;

addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
clc;
set(0, 'DefaultLineLineWidth', 1.5); set(0,'defaultAxesFontSize',16); 

rng default;
addpath('Functions');

%% Model parameters

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

Xu=-0.1068;

Xq=0.1192;

Mu=-5.9755;

Mq=-2.6478;

Xd=-10.1647;

Md=450.71;

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

% Noise
C_qax=[0, 1, 0; Xu, Xq, 0];
D_qax=[0; Xd];
sys_real = ss(A,B,C_qax,D_qax);

%%
%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                                %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.0076);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.01);                   %[rad/s]

% Delays
 
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

% Just used for time array
load ExcitationM

SetPoint=[0,0];

%% Values selected
simulation_time= ExcitationM(end,1);
 
%% Cutting the output

ta = 0: sample_time : simulation_time;

t = ta(1:end-1);

%% Imput

attenuation_sweep = 0.1;
start_freq = 0.2; % starting sweep frequency Hz
end_freq = 1.1; % ending sweep frequency Hz
start_time = 20; % [s] initial time of the sweep signal
end_time = 100; % [s] ending time of the sweep signal
phase = 0 ; % [rad] starting phase ot the input sweep 
[sw_custom_fcn] = variable_sweep(start_freq,end_freq,start_time,end_time,phase,t);
sw_custom = zeros(length(t),2);
sw_custom(:,1) = t;
sw_custom(:,2) = attenuation_sweep*sw_custom_fcn;
sw_pips = PIPS(sw_custom(:,2)); % compute the performance index of the sweep
fcn = sw_custom(:,2);
t_plot = sw_custom(:,1);

%% Simulation
simout = sim ('Simulator_Single_Axis');

%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

%% Estimation Bis

% Extract data from simulation
q = simout.q(1:end-1);
ax = simout.ax(1:end-1);
M = simout.Mtot(1:end-1);

u = M;

figure
subplot(411)
plot(t_plot, fcn)
hold on
grid
title('Excitation input')
subplot(412)
plot(t,M(1:length(t)))
hold on
grid
title('Normalised pitching moment') % From outside it's a sweep but the control distorts it
subplot(413)
plot(t,q(1:length(t)))
grid
title('Pitch rate [rad/s]')
subplot(414)
plot(t,ax(1:length(t)))
grid
title('Longitudinal acceleration [m/s^2]')
xlabel('Time [s]')

%% Iddata

data = iddata([q ax],u,sample_time);
data_fft = fft(data);

%% Model definition

odefun = @Model;

Xu = 0;
Xq = 0;
Mu = 0;
Mq = 0;
Xd = 0;
Md = 0;

parameters = {'Xu',Xu;'Xq',Xq;'Mu',Mu;'Mq',Mq;'X_delta',Xd;'M_delta', Md};
fcn_type = 'c';

sys = idgrey(odefun,parameters,fcn_type);

%% Grey-box estimation

opt = greyestOptions('Display','On');
sys = greyest(data_fft,sys,opt);
pvec = getpvec(sys);

cov_sys = getcov(sys);
diagonal = 3*sqrt(diag(cov_sys));
Err = trace(cov_sys);

%% 3-sigma uncertainties

Xu_u = ureal('Xu',pvec(1),'PlusMinus',[-diagonal(1) diagonal(1)]);
Xq_u = ureal('Xq',pvec(2),'PlusMinus',[-diagonal(2) diagonal(2)]);
Mu_u = ureal('Mu',pvec(3),'PlusMinus',[-diagonal(3) diagonal(3)]);
Mq_u = ureal('Mq',pvec(4),'PlusMinus',[-diagonal(4) diagonal(4)]);
Xd_u = ureal('Xd',pvec(5),'PlusMinus',[-diagonal(5) diagonal(5)]);
Md_u = ureal('Md',pvec(6),'PlusMinus',[-diagonal(6) diagonal(6)]);

g = 9.81;
A2 = [Xu_u Xq_u -g;
    Mu_u  Mq_u 0;
    0 1 0];
B2 = [Xd_u;
    Md_u;
    0];
C2 = [0 1 0;
     Xu_u Xq_u 0]; 

D2 = [0;Xd_u];

usys = ss(A2,B2,C2,D2);
figure
bode(usys,'b',usys.Nominal,'r--');
hold on 
bode(sys_real,'k:');
grid 
set(findall(gcf,'type','line'),'linewidth',2.5)
legend('Uncertainties','Estimated Parameters','Real System',Location="best")

%% Bode Error

[real_mag,real_phs,real_wout] = bode(sys_real);
real_Magnitude = squeeze(real_mag);
real_Phase = squeeze(real_phs);

[usys_mag,usys_phs,usys_wout] = bode(usys.Nominal);
usys_Magnitude = squeeze(usys_mag);
usys_Phase = squeeze(usys_phs);

figure
subplot(4,1,1)
semilogx(real_wout,abs(20*log10(real_Magnitude(1,:))-20*log10(usys_Magnitude(1,:))))
title('Error in Magnitude - q (dB)')
grid
subplot(4,1,2)
semilogx(real_wout,abs(real_Phase(1,:)-usys_Phase(1,:)))
title('Error in Phase - q (°)')
grid
subplot(4,1,3)
semilogx(real_wout,abs(20*log10(real_Magnitude(2,:))-20*log10(usys_Magnitude(2,:))))
title('Error in Magnitude - ax (dB)')
grid
subplot(4,1,4)
semilogx(real_wout,abs(real_Phase(2,:)-usys_Phase(2,:)))
title('Error in Phase - ax (°)')
grid
set(findall(gcf,'type','line'),'linewidth',1.5)

%% END OF CODE