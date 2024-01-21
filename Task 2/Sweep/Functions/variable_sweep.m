function [sweep] = variable_sweep(start_freq,end_freq,start_time,end_time,phase,t)
%VARIABLE_SWEEP generate a variable sweep with custom parameters
%   
% input
% start_freq        [rad/s] starting sweep frequency
% end_freq          [rad/s] ending sweep frequency
% start_time        [s] initial time of the sweep signal
% end_time          [s] ending time of the sweep signal
% phase             [rad] starting phase ot the input sweep
% t                 [s] time array where to print the sweep
%
% ouput
% sweep             custom sweep 

if start_time>end_time
    printf('Sweep function error:\n start time is bigger than end time');
end

custom_sweep = @(time) 1*sin(2*pi*(start_freq+(end_freq-start_freq).*(time-start_time)./(end_time-start_time)).*time+phase);
% sweep = custom_sweep(t);
sweep = custom_sweep(t).*(t>start_time).*(t<end_time);
end