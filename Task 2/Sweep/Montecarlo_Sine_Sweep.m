% To see the results go to "Results_plot" script.
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

%% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 

load ExcitationM

SetPoint=[0,0];

%% Values selected

simulation_time= ExcitationM(end,1);

%% Cutting the output

ta = 0: sample_time : simulation_time;

t = ta(1:end-1);

%% Diferrent inputs
% attenuation of the signal

% variable sweep input 
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

%% usample

n_MC = 100;

Xu_u = ureal('Xu',Xu,'percent',5);
Xq_u = ureal('Xq',Xq,'percent',5);
Mu_u = ureal('Mu',Mu,'percent',5);
Mq_u = ureal('Mq',Mq,'percent',5);
Xd_u = ureal('Xd',Xd,'percent',5);
Md_u = ureal('Md',Md,'percent',5);

Xu_sample = usample(Xu_u,n_MC);
Xq_sample = usample(Xq_u,n_MC);
Mu_sample = usample(Mu_u,n_MC);
Mq_sample = usample(Mq_u,n_MC);
Xd_sample = usample(Xd_u,n_MC);
Md_sample = usample(Md_u,n_MC);

Err_vec = zeros(n_MC,1);

param_estimated = cell(n_MC,1);
cov_param = cell(n_MC,1);
MONTECARLO_SW.err = cell(n_MC,1);
MONTECARLO_SW.err_rel = cell(n_MC,1);

%% MONTECARLO

for jj = 1: n_MC
    
    % Simulation
    Xu = Xu_sample(jj);
    Xq = Xq_sample(jj);
    Mu = Mu_sample(jj);
    Mq = Mq_sample(jj);
    Xd = Xd_sample(jj);
    Md = Md_sample(jj);
    
    MONTECARLO_SW.param{jj} = [Xu, Xq, Mu, Mq, Xd, Md]';  
    
    A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];
    
    B=[Xd; Md; 0];
    
    C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 
    
    D=[0; 0; 0; Xd];
    
    simout = sim ('Simulator_Single_Axis');
    
    % Delete temporary files
    
    if exist('slprj','dir')
        rmdir('slprj', 's')                                                    
    end
    
    % Estimation Bis
    
    % Extract data from simulation
    q = simout.q(1:end-1);
    ax = simout.ax(1:end-1);
    M = simout.Mtot(1:end-1);
    u = M;

    % Model definition
    odefun = @Model;

    Xu_ = 0;
    Xq_ = 0;
    Mu_ = 0;
    Mq_ = 0;
    Xd_ = 0;
    Md_ = 0;
    
    parameters = {'Xu',Xu_;'Xq',Xq_;'Mu',Mu_;'Mq',Mq_;'X_delta',Xd_;'M_delta', Md_};
    fcn_type = 'c';
    
    sys = idgrey(odefun,parameters,fcn_type);

    
    data = iddata([q ax],u,sample_time);
    data_fft = fft(data);
    
    % Grey-box estimation  
    opt = greyestOptions('Display','Off');
    sys = greyest(data_fft,sys,opt);
    pvec = getpvec(sys);
    cov_sys = getcov(sys);   
    
    MONTECARLO_SW.param_estimated{jj} = pvec;
    MONTECARLO_SW.cov_param{jj} = diag(cov_sys);    
    MONTECARLO_SW.trace(jj) = trace(cov_sys);
    
    MONTECARLO_SW.err{jj} = abs(pvec - MONTECARLO_SW.param{jj});
    MONTECARLO_SW.err_rel{jj} = MONTECARLO_SW.err{jj}./abs(MONTECARLO_SW.param{jj});  
end

%% Save
 
% save('MONTECARLO_SW','MONTECARLO_SW')

%% END OF CODE