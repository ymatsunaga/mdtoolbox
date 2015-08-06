function [lifetime, time_array] = calclifetime(data, recrossing_step)
%% calclifetime
% calculate lifetime by the maximum-likelihood estimation (MLE)
%
%% Syntax
%# lifetime = calcdlifetime(data)
%# lifetime = calcdlifetime(data, recrossing_step)
%# [lifetime, time_array] = calcdlifetime(data, recrossing_step)
%
%% Description
% lifetime of 1 in data is estimated by maximum likelihood
%
% * data - logical data true for the state of which the lifetime is calclated
% * recrossing_step - threshold step to ignore as recrossing events
%
%% Example
%#
% 
%% See also
%
%% References
% http://en.wikipedia.org/wiki/Poisson_distribution
% 

%% change point detection
% 1 for false->true, -1 for true-false, in 1 step evolution
index_change = data(2:end) - data(1:end-1);

% just before inward crossing (false->true)
index_inward = find(index_change == 1);

% just before outward crossing (true->false)
index_outward = find(index_change == -1);

%% extrude exceptional cases
% exclude the case where the initial step is already in true
if index_outward(1) < index_inward(1)
  index_outward(1) = [];
end

% exclude the case where the final step remains in true
if index_outward(end) < index_inward(end)
  index_inward(end) = [];
end

assert(numel(index_inward) == numel(index_outward), ['inward and outward crossing events do not match']);

%% calculate crossing time
ncrossing = numel(index_inward);
time_array = zeros(ncrossing, 1);

for i = 1:ncrossing
  time_array(i) = index_outward(i) - index_inward(i) - 1;
end

if nargin > 1
  lifetime = mean(time_array(time_array > recrossing_step));
else
  lifetime = mean(time_array);
end

