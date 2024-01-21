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
addpath('Save')
clc;
set(0, 'DefaultLineLineWidth', 1.5); set(0,'defaultAxesFontSize',16); 

rng default;

%% Selection

% Run different Randomic Binary Sequence ([ON,OFF])
run_RBS = 'OFF';
% Number of RBS generated (If run_RBS = 'ON')
N_RBS = 50;

% Montecarlo run to test sensitivity of the chosen RBS ([ON,OFF])
run_MC = 'OFF';
% Number of MC samples considered (If run_MC = 'ON')
n_MC = 100;

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

% noise.Enabler = 0;
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

%% Real parameters
preal = [-0.1068;
          0.1192;
         -5.9755;
         -2.6478;
        -10.1647;
        450.71];

%% Model definition

odefun = @Model;

% Initial guess for the parameters
Xu = 0;
Xq = 0;
Mu = 0;
Mq = 0;
Xd = 0;
Md = 0;

Err_ExcitationM = 0.000871854135723676; % As computed in task 1

parameters = {'Xu',Xu;'Xq',Xq;'Mu',Mu;'Mq',Mq;'X_delta',Xd;'M_delta', Md};
fcn_type = 'c';

sys_0 = idgrey(odefun,parameters,fcn_type);

%% RBS input

col =  [0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.4660 0.6740 0.1880;
        0.6350 0.0780 0.1840];

t = ExcitationM(:,1);
simulation_time = t(end)-t(1);

% Random Binary Sequence
attenuation = 0.2;
t_0 = 0; % starting time of the RBS
duration = 60; %[s] duration time of the signal
T = 0.25; % time of switching states / sampling interval
t_rbs = 0 : sample_time : (t_0 + duration + 15);
t = t_rbs;

if strcmp(run_RBS,'ON')
    
    RBS_input = cell(N_RBS,1);
    PIPS_value = zeros(N_RBS,1);
    est_err = zeros(6,N_RBS);
    est_err_rel = zeros(6,N_RBS);
    Err = zeros(N_RBS,1);
    kk = 1;
    best_RBS_error = 1;
    best_kk = 1;

    f = waitbar(0, 'Starting');        
    for jj = 1:N_RBS
    
        waitbar(jj/N_RBS, f, sprintf('Progress: %d %%', floor(jj/N_RBS*100)));       
    
        [y_RBS_fcn,t_rbs] = dis2time_RBS(t_0,duration,T,t_rbs); % discrete
        y_RBS = zeros(length(t_rbs),2);
        y_RBS(:,1) = t_rbs;
        y_RBS(:,2) = attenuation*y_RBS_fcn;
        
        RBS_input{jj}.RBS = y_RBS;
        
        PIPS_value(jj) = PIPS(y_RBS(:,2)); % compute the performance index of the RBS
        RBS_input{jj}.PIPS = PIPS_value(jj);
        
        % Simulation
        simout = sim ('Simulator_Single_Axis_RBS');
        
        % Extract data from simulation
        q = simout.q(1:end);
        ax = simout.ax(1:end);
        M = simout.Mtot(1:end);

        RBS_input{jj}.q = q;
        RBS_input{jj}.ax = ax;
        RBS_input{jj}.M = M;
        
        u = M;
        
        data = iddata([q ax],u,sample_time);
        data_fft = fft(data);
        
        % Grey-box estimation
        sys = greyest(data_fft,sys_0);
        
        pvec = getpvec(sys);
        RBS_input{jj}.pvec = pvec;
        RBS_input{jj}.cov_sys = getcov(sys);
    
        % Absolute estimation errors
        est_err(:,jj) = abs(preal - pvec);
        % Relative estimation errors
        est_err_rel(:,jj) = est_err(jj)./abs(preal);
    
        RBS_input{jj}.est_err = est_err(:,jj);
        RBS_input{jj}.est_err_rel = est_err_rel(:,jj);
    
        Err(jj) = trace(RBS_input{jj}.cov_sys);
        RBS_input{jj}.Err = Err(jj);
    
        if Err_ExcitationM > Err(jj)
            best_sequence{kk}.err = Err(jj);
            best_sequence{kk}.best_RBS = y_RBS;
            best_sequence{kk}.pvec = pvec;
            best_sequence{kk}.cov_sys = getcov(sys);
            best_sequence{kk}.jj = jj;
           
            if best_sequence{kk}.err < best_RBS_error
                best_RBS_error = best_sequence{kk}.err;
                best_kk = kk;
            end

            kk = kk+1;
        end
        pause(0.01)  
    end
    close(f)
    
    save('best_RBS','best_sequence','RBS_input','PIPS_value','Err','best_kk')
    
    % Best RBS local
    best_sequence{best_kk}.err
    best_sequence{best_kk}.pvec
    best_sequence{best_kk}.cov_sys

    % Mean PIPS and variance
    mean_PIPS = mean(PIPS_value);
    var_PIPS = var(PIPS_value);
    
    % Mean error and variance
    mean_err = mean(Err);
    var_err = var(Err);
    
    % Delete temporary files
    if exist('slprj','dir')
        rmdir('slprj', 's')                                                    
    end
    
    else
    load('best_RBS')
end

%% Plot traces of covariance matrices vs trace error ExcitationM

figure
set(gca, 'YScale', 'log')
hold on
grid minor
plot(1:N_RBS,Err_ExcitationM*ones(N_RBS,1),'Marker','.','MarkerSize',15,'Color',col(2,:),'LineWidth',2,'DisplayName','ExcitationM trace value')
scatter(1,RBS_input{1}.Err,15,'MarkerEdgeColor',col(4,:),'MarkerFaceColor',col(4,:),'LineWidth',1.5,'DisplayName','RBS trace values')

for jj = 2:N_RBS
    scatter(jj,RBS_input{jj}.Err,15,'MarkerEdgeColor',col(4,:),'MarkerFaceColor',col(4,:),'LineWidth',1.5,'HandleVisibility','off')
end

xlabel('Sample','Interpreter','latex')
ylabel('Trace of covariance matrix','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');
legend('Location','best','Interpreter','latex')

%% Sensitivity of the input sequence

if strcmp(run_MC,'ON')

    % Real parameters
    Xu = -0.1068;
    Xq = 0.1192;
    Mu = -5.9755;
    Mq = -2.6478;
    Xd = -10.1647;
    Md = 450.71;
    
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
    
    n_RBS = length(best_sequence);
    RBS_MC = cell(n_MC,n_RBS);
    MC_out = cell(size(best_sequence));
    est_err = zeros(6,n_MC);
    est_err_rel = zeros(6,n_MC);
    Err = zeros(n_MC,1);
            
    for kk = 1:n_RBS
    
        y_RBS = best_sequence{kk}.best_RBS;
    
        for jj = 1: n_MC
            
            X = ['Iteration Number ',num2str(jj),' of ',num2str(n_MC),' for the sequence ',num2str(kk),' of ',num2str(n_RBS)];
            disp(X)
        
            % Simulation
            Xu = Xu_sample(jj);
            Xq = Xq_sample(jj);
            Mu = Mu_sample(jj);
            Mq = Mq_sample(jj);
            Xd = Xd_sample(jj);
            Md = Md_sample(jj);
    
            preal_MC = [Xu;Xq;Mu;Mq;Xd;Md];
            RBS_MC{jj,kk}.preal = preal_MC;
        
            A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];
            B=[Xd; Md; 0];
            C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0];
            D=[0; 0; 0; Xd];
        
            % Simulation
            simout = sim ('Simulator_Single_Axis_RBS');
        
            % Extract data from simulation
            q = simout.q(1:end);
            ax = simout.ax(1:end);
            M = simout.Mtot(1:end);
            
            u = M;
            
            data = iddata([q ax],u,sample_time);
            data_fft = fft(data);
        
            % Grey-box estimation
            sys = greyest(data_fft,sys_0);
        
            pvec = getpvec(sys);
            RBS_MC{jj,kk}.pvec = pvec;        
            RBS_MC{jj,kk}.cov_sys = getcov(sys);
        
            % Absolute estimation errors
            est_err(:,jj) = abs(preal_MC - pvec);
            % Relative estimation errors
            est_err_rel(:,jj) = abs(est_err(jj)./preal_MC);
        
            RBS_MC{jj,kk}.est_err = est_err(:,jj);
            RBS_MC{jj,kk}.est_err_rel = est_err_rel(:,jj);
        
            Err(jj) = trace(RBS_MC{jj}.cov_sys);
            RBS_MC{jj,kk}.Err = Err(jj);
        end
    
        MC_out{kk}.Err_mean = mean(Err);
        MC_out{kk}.Err_var = var(Err);
        MC_out{kk}.Err_pvec = mean(abs(est_err),2);
        MC_out{kk}.Err_rel_pvec = mean(abs(est_err_rel),2);

        if kk ~= 1
            if best_est > MC_out{kk}.Err_mean
                best_kk_MC = kk;
                best_est = MC_out{kk}.Err_mean;
            end
        else
            best_kk_MC = 1;
            best_est = MC_out{kk}.Err_mean;
        end
    end
    
    save('bestMC','RBS_MC','MC_out','best_kk_MC')
    
    % Delete temporary files
    if exist('slprj','dir')
        rmdir('slprj', 's')                                                    
    end

else
    load('bestMC')
end

%% Bode plot

pvec = best_sequence{best_kk_MC}.pvec;
dcov = 3*sqrt(diag(best_sequence{best_kk_MC}.cov_sys));

Xu_u = ureal('Xu',pvec(1),'PlusMinus',[-dcov(1),dcov(1)]);
Xq_u = ureal('Xq',pvec(2),'PlusMinus',[-dcov(2),dcov(2)]);
Mu_u = ureal('Mu',pvec(3),'PlusMinus',[-dcov(3),dcov(3)]);
Mq_u = ureal('Mq',pvec(4),'PlusMinus',[-dcov(4),dcov(4)]);
Xd_u = ureal('Xd',pvec(5),'PlusMinus',[-dcov(5),dcov(5)]);
Md_u = ureal('Md',pvec(6),'PlusMinus',[-dcov(6),dcov(6)]);

Au = [Xu_u, Xq_u, -9.81; Mu_u, Mq_u, 0; 0, 1, 0];
Bu = [Xd_u; Md_u; 0];
Cu = [0, 1, 0; Xu_u, Xq_u, 0];
Du = [0; Xd_u];

figure
usys = ss(Au,Bu,Cu,Du);
bode(usys,'b',usys.Nominal,'r--');
legend(Location="best")
grid
set(findall(gcf,'type','line'),'linewidth',2.5)

% Real system
Xu = -0.1068;
Xq = 0.1192;
Xd = -10.1647;

C_qax = [0, 1, 0; Xu, Xq, 0];
D_qax = [0; Xd];

figure
sys_real = ss(A,B,C_qax,D_qax);
bode(usys,'b',usys.Nominal,'r--');
hold on
bode(sys_real,'k:');
grid
set(findall(gcf,'type','line'),'linewidth',2.5)
legend('Uncertainties','Estimated Parameters','Real System',Location="best")

%% Bode plot Error

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

%% RBS with the best sensitivity

figure
subplot(411)
plot(best_sequence{best_kk_MC}.best_RBS(:,1),best_sequence{best_kk_MC}.best_RBS(:,2),'LineWidth',2)
axis on
grid on
ylim([-0.05 0.25])
xlim([0 75])
title('Excitation Input','FontSize',16)
subplot(412)
plot(0:4e-3:simulation_time,RBS_input{best_sequence{best_kk_MC}.jj}.M)
hold on
grid
xlim([0 75])
title('Normalised pitching moment') % From outside it's a sweep but the control distorts it
subplot(413)
plot(0:4e-3:simulation_time,RBS_input{best_sequence{best_kk_MC}.jj}.q)
grid
xlim([0 75])
title('Pitch rate [rad/s]')
subplot(414)
plot(0:4e-3:simulation_time,RBS_input{best_sequence{best_kk_MC}.jj}.ax)
grid
xlim([0 75])
title('Longitudinal acceleration [m/s^2]')
xlabel('Time [s]')

%% RBS relative errors

figure
set(gca, 'YScale', 'log')
hold on
grid minor

preal = [-0.1068;
          0.1192;
         -5.9755;
         -2.6478;
        -10.1647;
        450.71];

est_rel = abs(best_sequence{best_kk_MC}.pvec-preal);

for rr = 1:6
    plot(1:n_MC,est_rel(rr)*ones(n_MC,1),'Marker','.','MarkerSize',15,'Color',col(rr,:),'LineWidth',2)
end

for jj = 1:n_MC
    for rr = 1:6
        scatter(jj,abs(RBS_MC{jj,best_kk_MC}.pvec(rr)-RBS_MC{jj,best_kk_MC}.preal(rr)),15,'MarkerEdgeColor',col(rr,:),'MarkerFaceColor',col(rr,:),'LineWidth',1.5,'HandleVisibility','off')
    end
end

legend('Xu','Xq','Mu','Mq','Xd','Md','Location','bestoutside','Interpreter','latex')
xlabel('Sample','Interpreter','latex')
ylabel('Error','Interpreter','latex')
title('Relative error on each parameter','FontSize',16)
set(gca,'TickLabelInterpreter','latex');

%% Check on the variation of the parameters

figure
set(gca, 'YScale', 'log')
hold on
grid minor
for jj = 1:n_MC
    variation = abs(RBS_MC{jj,best_kk_MC}.preal-preal)./abs(preal);
    for rr = 1:6
        scatter(jj,variation(rr),15,'MarkerEdgeColor',col(rr,:),'MarkerFaceColor',col(rr,:),'LineWidth',1.5,'HandleVisibility','off')
    end
end

xlabel('Sample','Interpreter','latex')
ylabel('Variation of the real parameters','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');

%% Plot trace of covariance matrix

figure
set(gca, 'YScale', 'log')
hold on
grid minor
for jj = 1:n_MC
    scatter(jj,RBS_MC{jj,best_kk_MC}.Err,15,'MarkerEdgeColor',col(4,:),'MarkerFaceColor',col(4,:),'LineWidth',1.5,'HandleVisibility','off')
end

xlabel('Sample','Interpreter','latex')
ylabel('Trace of covariance matrix','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');

%% END OF CODE