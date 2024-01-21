function [Performance_index] = PIPS(signal)
%PIPS return a key performance index of a perturbation signal
%   
% input
%   signal              signal to be analized
%
% output
%   performance index   the resulting PIPS with a relative score from
%                       0-100 %

u_rms = rms(signal); % root mean square of the signal
u_mean = mean(signal); % mean of the signal
u_max = max(signal);
u_min = min(signal);

Performance_index = 200*((u_rms^2-u_mean^2)^0.5)/(u_max-u_min);
end