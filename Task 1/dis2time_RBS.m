function [y_RBS,t] = dis2time_RBS(t_0,duration,T,t)
%DIS2TIME Generate a Random Binary Sequence in time domain
%
% input
%   t_0         [s] starting time of the RBS
%   duration    [s] duration time of the signal
%   T           [s] time of switching states / sampling interval
%   t           [s] time array where to print the RBS (uniform time array)
%
% output
%   y_RBS       Random Binary Sequence in time domain
%   t           [s] time array of the RBS

t_RBS = 0;% initialising time of the RBS
s1 = round(rand(1,1)); % initialising of the state variable
N = 5; % integer of rounding
sample_time = round(t(end)/length(t),N);
y_RBS = zeros(1,length(t)); 
t_real = t_0+t_RBS; % time in the simulation
index1 = floor(t_real/sample_time)+1; % index in the time array
index2 = floor((t_real+T)/sample_time); % index of the end interval in the
                                    % time array
while (t_RBS < duration)
    y_RBS(1,index1:index2) = s1 ;
    t_RBS = t_RBS+T; % time of the switchinf state
    t_real = t_0 + t_RBS; % time ragrding the simulation
    index1 = index2+1;
    index2 = floor((t_real+T)/sample_time);
    s1 = round(rand(1,1)); % random number 0/1 generator
end

end