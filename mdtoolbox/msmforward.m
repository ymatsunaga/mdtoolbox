function [logL, alpha, factor] = msmforward(data, T, emission, pi_i)
%% msmforward
% forward algorithm to calculate the probability of observed sequence
%
%% Syntax
%# [logL, alpha] = msmforward(data, T, emission, pi_i)
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
if iscell(data)
  data_cell = data;
  is_cell   = true;
else
  data_cell{1} = data;
  is_cell      = false;
end

ndata = numel(data_cell);
nstate = size(T, 1);

if ~exist('emission', 'var') || isempty(emission)
  emission = eye(nstate);
end
emission = emission';

if iscolumn(pi_i)
  pi_i = pi_i';
end

%% forward algorithm
logL        = zeros(ndata, 1);
alpha_cell  = cell(ndata, 1);
factor_cell = cell(ndata, 1);

for idata = 1:ndata
  data = data_cell{idata};
  nframe = numel(data);
  alpha  = zeros(nframe, nstate);
  factor = zeros(nframe, 1);
  
  alpha(1, :) = pi_i.*emission(data(1), :);
  factor(1) = sum(alpha(1, :));
  alpha(1, :) = alpha(1, :)./factor(1);

  for iframe = 2:nframe
    alpha(iframe, :) = sum(bsxfun(@times, alpha(iframe-1, :)', T)) .* emission(data(iframe), :);
    factor(iframe) = sum(alpha(iframe, :));
    alpha(iframe, :) = alpha(iframe, :)./factor(iframe);
  end

  logL(idata) = sum(log(factor));
  alpha_cell{idata} = alpha;
  factor_cell{idata} = factor;
end

if is_cell
  alpha  = alpha_cell;
  factor = factor_cell;
else
  alpha  = alpha_cell{1};
  factor = factor_cell{1};
end

