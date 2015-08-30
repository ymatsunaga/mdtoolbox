function [logL, beta, beta_nonscaled] = msmbackward(data, factor, T, emission, pi_i)
%% msmbackward
% backward algorithm to calculate the probability of sequence of observarion
%
%% Syntax
%# [logL, beta, beta_nonscaled] = msmbackward(data, factor, T, emission, pi_i)
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
nframe = numel(data);
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
  beta(iframe, :) = sum(bsxfun(@times, T, emission(:, data(iframe+1)) .* beta(iframe+1, :)'))./factor(iframe+1);
end

logL = sum(log(factor));

if nargout > 2
  beta_nonscaled = beta;
  for iframe = 1:(nframe-1)
    beta_nonscaled(iframe, :) = prod(factor(iframe+1:nframe))*beta(iframe, :);
  end
end

