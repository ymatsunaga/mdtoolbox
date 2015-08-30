function [logL, alpha, factor] = msmforward(data, T, emission, pi_i)
%% msmforward
% forward algorithm to calculate the probability of sequence of observarion
%
%% Syntax
%# [logL, alpha, factor] = msmforward(data, T, emission, pi_i)
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
emission = emission';

alpha = zeros(nframe, nstate);
factor = zeros(nframe, 1);

%% forward algorithm
alpha(1, :) = pi_i.*emission(data(1), :);
factor(1) = sum(alpha(1, :));
alpha(1, :) = alpha(1, :)./factor(1);

for iframe = 2:nframe
  alpha(iframe, :) = sum(bsxfun(@times, alpha(iframe-1, :)', T)) .* emission(data(iframe), :);
  factor(iframe) = sum(alpha(iframe, :));
  alpha(iframe, :) = alpha(iframe, :)./factor(iframe);
end

logL = sum(log(factor));

