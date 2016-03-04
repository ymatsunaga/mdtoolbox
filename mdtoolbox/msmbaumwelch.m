function [T, emission, pi_i] = msmbaumwelch(data, T0, emission0, pi0_i)
%% msmbaumwelch
% Baum-Welch algorithm
%
%% Syntax
%# [T, emission, pi_i] = msmbaumwelch(data, T0, emission0, pi0_i)
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

%% setup
if iscell(data)
  data_cell = data;
else
  data_cell{1} = data;
end

%nframe = numel(data);
ndata = numel(data_cell);
nstate = size(T0, 1);
nobs = size(emission0, 2);

if iscolumn(pi0_i)
  pi0_i = pi0_i';
end

TOLERANCE = 10^(-4);
check_convergence = inf;
count_iteration = 0;
logL_old = 1.0;
while check_convergence > TOLERANCE
  %% E-step
  [logL, alpha, factor] = msmforward(data_cell, T0, emission0, pi0_i);
  [~, beta] = msmbackward(data_cell, factor, T0, emission0, pi0_i);

  log_alpha     = cellfun(@log, alpha, 'UniformOutput', false);
  log_beta      = cellfun(@log, beta,  'UniformOutput', false);
  log_T0        = log(T0);
  log_emission0 = log(emission0);

  %% M-step
  % equilibrium probability (pi_i)
  pi_i = zeros(1, nstate);
  log_gamma = cell(ndata, 1);
  for idata = 1:ndata
    log_gamma{idata} = log_alpha{idata} + log_beta{idata};
    pi_i = pi_i + exp(log_gamma{idata}(1, :));
  end
  pi_i  = pi_i./sum(pi_i);

  % emission probability (emission)
  emission = zeros(nstate, nobs);
  for idata = 1:ndata
    data = data_cell{idata};
    for istate = 1:nstate
      for iobs = 1:nobs
        id = find(data == iobs);
        emission(istate, iobs) = emission(istate, iobs) + sum(exp(log_gamma{idata}(id, istate)));
      end
    end
  end
  emission = bsxfun(@rdivide, emission, sum(emission, 2));
  emission(isnan(emission)) = 0;

  % transition probability (T)
  T = zeros(nstate, nstate);
  for idata = 1:ndata
    data = data_cell{idata};
    nframe = numel(data);
    for iframe = 2:nframe
      log_xi = bsxfun(@plus, log_alpha{idata}(iframe-1, :)', log_beta{idata}(iframe, :));
      T = T + exp(bsxfun(@plus, log_xi, log_emission0(:, data(iframe))') + log_T0)./factor{idata}(iframe);
    end
  end
  T = bsxfun(@rdivide, T, sum(T, 2));
  T(isnan(T)) = 0;

  % check convergence
  count_iteration = count_iteration + 1;
  logL = sum(logL);
  check_convergence = abs(logL_old - logL);
  if mod(count_iteration, 1) == 0
    fprintf('%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n', count_iteration, logL, check_convergence, TOLERANCE);
  end

  logL_old = logL;
  pi0_i = pi_i;
  emission0 = emission;
  T0 = T;
end

