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
nframe = numel(data);
nstate = size(T0, 1);
nobs = size(emission0, 2);

if iscolumn(pi0_i)
  pi0_i = pi0_i';
end

TOLERANCE = 10^(-8);
check_convergence = inf;
count_iteration = 0;
logL_old = 1.0;
while check_convergence > TOLERANCE
  %% E-step
  [logL, alpha, factor] = msmforward(data, T0, emission0, pi0_i);
  [~, beta] = msmbackward(data, factor, T0, emission0, pi0_i);

  log_alpha     = log(alpha);
  log_beta      = log(beta);
  log_T0        = log(T0);
  log_emission0 = log(emission0);

  %% M-step
  % pi_i
  log_gamma = log_alpha + log_beta;
  pi_i = exp(log_gamma(1, :));
  pi_i  = pi_i./sum(pi_i);

  % emission
  emission = zeros(nstate, nobs);
  for istate = 1:nstate
    for iobs = 1:nobs
      id = find(data == iobs);
      emission(istate, iobs) = emission(istate, iobs) + sum(exp(log_gamma(id, istate)));
    end
  end
  emission = bsxfun(@rdivide, emission, sum(emission, 2));
  emission(isnan(emission)) = 0;

  % T
  T = zeros(nstate, nstate);
  for iframe = 2:nframe
    log_xi = bsxfun(@plus, log_alpha(iframe-1, :)', log_beta(iframe, :));
    T = T + exp(bsxfun(@plus, log_xi, log_emission0(:, data(iframe))') + log_T0)./factor(iframe);
  end
  T = bsxfun(@rdivide, T, sum(T, 2));
  T(isnan(T)) = 0;

  % check convergence
  count_iteration = count_iteration + 1;
  abs_logL = abs(logL);
  check_convergence = abs(logL_old - logL);
  if mod(count_iteration, 1) == 0
    fprintf('%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n', count_iteration, logL, check_convergence, TOLERANCE);
  end

  logL_old = logL;
  pi0_i = pi_i;
  emission0 = emission;
  T0 = T;
end

