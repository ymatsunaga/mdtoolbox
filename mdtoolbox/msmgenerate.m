function [state, observation] = msmgenerate(nframe, T, emission, pi_i)
%% msmgenerate
% generate discrete Markov state trajectory
%
%% Syntax
%# [state, observation] = msmgenerate(nframe, T, emission, pi_i)
%
%% Description
% 
%
%% Example
%#
% 
%% See also
%
%% TODO
% 
%

state = zeros(nframe, 1);
observation = zeros(nframe, 1);

state(1)       = sample(pi_i);
%state(1)
observation(1) = sample(emission(state(1), :));
%observation(1)
for iframe = 2:nframe
  %iframe
  state(iframe)       = sample(T(state(iframe-1), :));
  observation(iframe) = sample(emission(state(iframe), :));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = sample(p)
%disp('called')
i = min(find(rand <= cumsum(p) / sum(p)));
%disp('called2')
%numel(i)
if numel(i) > 1
  i = i(1);
end

