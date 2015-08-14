function [pmf_m, t_k] = dtram(c_k, fhandle_k, temperature, edge_m)
%% dtram
% calculate (reduced) potential of mean force and transition matrix by using dTRAM
%
%% Syntax
%# [pmf_m, t_k] = dtram(c_k, fhandle_k, temperature, edge_m)

%% Description
%
% * c_k         - transition count matrix for discrete Markov states
%                 [cell K x 1] [each cel has M x M matrix]
% * fhandle_k   - cell of function handles which represent bias-factors
%                 [cell K x 1]
% * temperature - Temperature in Kelvin
%                 [double scalar]
% * edge_m      - edges of state-bins
%                 [double 1 x M or M x 1]
%
% * pmf_m       - potential of mean force in data-bins under unbiased condition
%                 [double 1 x M]
% * t_k         - estimated transition matrix
%                 [double K x 1] [each cel has M x M matrix]
% 
%% Example
%# see http://mdtoolbox.readthedocs.org/en/latest/alat_1D_umbrella_wham.html
% 
%% See also
% wham msmcountmatrix
%
%% References
% This function uses the method described in
% [1]	H. Wu, A. S. J. S. Mey, E. Rosta, and F. Noé, 
%     J. Chem. Phys. 141, 214106 (2014).
%
%% TODO
% sparse matrix support
% mex for kernel part
%

%% preparation
% Boltzmann constant in kcal/(mol K)
KB = 0.00198719168260038;
% Tolerance for the convergence of iterations
TOLERANCE = 10^(-2);
%TOLERANCE = 10^(-2);

% K: number of thermodynamic states
K = numel(c_k);

% check consistency of the number of thermodynamic states
if K ~= numel(fhandle_k)
  error('# of thermodynamic states do not match... data has %d windows. fhandle has %d windows.', numel(data_k), numel(fhandle_k));
end

% M: number of Markov states
M = size(c_k{1}, 1);

% convert from sparse to full matrix
for k = 1:K
  if issparse(c_k{1})
    c_k{k} = full(c_k{k});
  end
end

% centers of bins
if exist('edge_m', 'var') && ~isempty(edge_m)
  center_m = 0.5 * (edge_m(2:end) + edge_m(1:(end-1)));
else
  center_m = [];
end

%% calculate bias-factor
% bias_km: bias-factor for umbrella-window (k) at bin-center (m)
bias_km = zeros(K, M);
for k = 1:K
  for m = 1:M
    bias_km(k, m) = fhandle_k{k}(center_m(m));
  end
end
% convert to dimensionless energies
bias_km = bias_km ./ (KB*temperature);

%% flat prior
c_k_cell = c_k;
clear c_k;
c_k = zeros(M, M, K);
for k = 1:K
  c_k(:, :, k) = c_k_cell{k} + (1./M)*10^(-5);
end

%% prepare symmetrized transition count matrix
c_sym_k = zeros(M, M, K);
for k = 1:K
  c_sym_k(:, :, k) = c_k(:, :, k) + c_k(:, :, k)';
end

%% neighbor prior
%% (in principle, MATLAB's log() exp() allow negative c_ij except for zero)
% for k = 1:K
%   ind_zero = find(c_sym_k{k}==0);
%   c_k{k}(ind_zero) = - 1;
% end
%
% c_sym_k = {};
% for k = 1:K
%   c_sym_k{k} = c_k{k} + c_k{k}';
% end

%% solve the dTRAM equations by self-consistent iteration
prior = double(eps('single'));
%prior = -1;
%prior = 0;

log_lambda_km = zeros(K, M); % log of Langrange's multipliers lambda
for k = 1:K
  log_lambda_km(k, :) = log(sum(c_k(:,:,k), 2)' + (M*prior) + 10^(-10));
  %log_lambda_km(k, :) = log(sum(c_k{k}, 2)' + 10^(-10));
end
check_convergence = inf;
t_k = log_lambda_km;

prob_m = ones(1, M)./M;
log_prob_m = log(prob_m);
pmf_m = - log_prob_m;
pmf_new_m = pmf_m;

log_lambda_new_km = log_lambda_km;
log_prob_new_m   = log_prob_m;

count_iteration = 0;
while check_convergence > TOLERANCE

  % 1st equation (Eq. 16 of Ref. 1)
  for k = 1:K
    log_iij_mn = bsxfun(@plus, log_lambda_km(k, :), (-bias_km(k,:)+log_prob_m)');
    log_iij_mn2 = real(log_iij_mn);
    log_jji_mn = bsxfun(@plus, (-bias_km(k,:)+log_prob_m), log_lambda_km(k, :)');
    log_jji_mn2 = real(log_jji_mn);
    vmax = max(log_iij_mn2, log_jji_mn2);
    log_iij_jji_mn = vmax + log(exp(log_iij_mn - vmax) + exp(log_jji_mn - vmax));

    log_c_mn = log(c_sym_k(:, :, k));
    log_mn = bsxfun(@plus, bsxfun(@plus, log_c_mn-log_iij_jji_mn, (-bias_km(k,:)+log_prob_m)), log_lambda_km(k,:)');
    log_lambda_new_km(k, :) = logsumexp2(log_mn');
  end
  log_lambda_km = log_lambda_new_km;

  % 2nd equation (Eq. 17 of Ref. 1)
  for i = 1:M
    log_numerator_km = log(squeeze(c_k(:, i, :))');
    
    log_iij_km = bsxfun(@plus, log_lambda_km, -bias_km(:, i)) + log_prob_m(i);
    log_iij_km2 = real(log_iij_km);
    log_jji_km = bsxfun(@plus, bsxfun(@plus, -bias_km, log_lambda_km(:, i)), log_prob_m);
    log_jji_km2 = real(log_jji_km);
    vmax = max(log_iij_km2, log_jji_km2);
    log_iij_jji_km = vmax + log(exp(log_iij_km - vmax) + exp(log_jji_km - vmax));

    log_c_km = log(squeeze(c_sym_k(i, :, :))');
    log_denominator_km = bsxfun(@plus, log_c_km + log_lambda_km - log_iij_jji_km, -bias_km(:, i));
    log_prob_new_m(i) = real(logsumexp(log_numerator_km(:)) - logsumexp(log_denominator_km(:)));
  end
  
  prob_new_m = exp(log_prob_new_m);
  prob_new_m = prob_new_m./sum(prob_new_m);
  pmf_new_m  = -log(prob_new_m);
  
  % check convergence
  check_convergence = max(abs(pmf_new_m - pmf_m))./std(pmf_new_m);

  prob_m = prob_new_m;
  log_prob_m = log(prob_new_m);
  pmf_m = pmf_new_m;

  count_iteration = count_iteration + 1;
  if mod(count_iteration, 10) == 0
    fprintf('%dth iteration  delta = %e  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
    fprintf('potential of mean force = ');
    fprintf('%f ', pmf_m);
    fprintf('\n');
    fprintf('\n');
  end
end

%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
% max_x = max(x);
% exp_x = exp(bsxfun(@minus, x, max_x));
% s = log(sum(exp_x)) + max_x;
max_x = max(real(x));
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;
%s = real(s);

function s = logsumexp2(x)
% max_x = max(x);
% exp_x = exp(bsxfun(@minus, x, max_x));
% s = log(sum(exp_x)) + max_x;
max_x = max(real(x));
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

