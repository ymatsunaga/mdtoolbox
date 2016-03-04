function [logL, beta, beta_nonscaled] = msmbackward(data, factor, T, emission, pi_i)
%% msmbackward
% backward algorithm to calculate the probability of observed sequence
%
%% Syntax
%# [logL, beta] = msmbackward(data, factor, T, emission, pi_i)
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
if iscell(data)
  data_cell   = data;
  factor_cell = factor;
  is_cell     = true;
else
  data_cell{1}   = data;
  factor_cell{1} = factor;
  is_cell        = false;
end

ndata = numel(data_cell);
nstate = size(T, 1);

if ~exist('emission', 'var') || isempty(emission)
  emission = eye(nstate);
end

if iscolumn(pi_i)
  pi_i = pi_i';
end

T = T';

%% backward algorithm
logL      = zeros(ndata, 1);
beta_cell = cell(ndata, 1);
if nargout > 2
  beta_nonscaled_cell = cell(ndata, 1);
end

for idata = 1:ndata
  data   = data_cell{idata};
  factor = factor_cell{idata};
  nframe = numel(data);
  beta   = zeros(nframe, nstate);

  beta(nframe, :) = 1;

  for iframe = (nframe-1):-1:1
    beta(iframe, :) = sum(bsxfun(@times, T, emission(:, data(iframe+1)) .* beta(iframe+1, :)'))./factor(iframe+1);
  end

  logL(idata) = sum(log(factor));
  beta_cell{idata} = beta;
  if nargout > 2
    beta_nonscaled = beta;
    for iframe = 1:(nframe-1)
      beta_nonscaled(iframe, :) = prod(factor(iframe+1:nframe))*beta(iframe, :);
    end
    beta_nonscaled_cell{idata} = beta_nonscaled;
  end
end

if is_cell
  beta = beta_cell;
  if nargout > 2
    beta_nonscaled = beta_nonscaled_cell;
  end
else
  beta = beta_cell{1};
  if nargout > 2
    beta_nonscaled = beta_nonscaled_cell{1};
  end
end

