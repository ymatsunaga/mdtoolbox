function [likelihood, beta, beta_nonscaled] = msmbackward(T, pi_i, emission, observation, factor)
%% msmbackward
% backward algorithm to calculate the probability of sequence of observarion
%
%% Syntax
%# [likelihood, beta] = msmbackward(T, pi_i, emission, observation, factor)
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
% robustness against zero elements of T
%

%% setup
nframe = numel(observation);
nstate = size(T, 1);

if iscolumn(pi_i)
  pi_i = pi_i';
end

if ~exist('emission', 'var') || isempty(emission)
  emission = eye(nstate);
end

T = T';

beta = zeros(nframe, nstate);

%% backward algorithm
beta(nframe, :) = 1;

for iframe = (nframe-1):-1:1
  beta(iframe, :) = sum(bsxfun(@times, T, emission(:, observation(iframe+1)) .* beta(iframe+1, :)'))./factor(iframe+1);
end

likelihood = prod(factor);

if nargout > 2
  beta_nonscaled = beta;
  for iframe = 1:(nframe-1)
    beta_nonscaled(iframe, :) = prod(factor(iframe+1:nframe))*beta(iframe, :);
  end
end

