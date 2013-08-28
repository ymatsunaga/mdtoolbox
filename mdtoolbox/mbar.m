function f_k = mbar(u_kn, fhandle_k, data_kn, tolerance)
%% mbar
% calculate the free energy differences of umbrella-windowed systems by using the Multistate Bennet Acceptance Ratio Method (MBAR)
%
%% Syntax
%# f_k = mbar(u_kn, fhandle_k, data_kn)
%# f_k = mbar(u_kn, fhandle_k, data_kn, tolerance)
%
%% Description
%
% * u_kn      - unbiased (dimensionless) potential energy of n-th snapshot from k-th umbrella-windows
%               [cell nwindow x 1]
% * fhandle_k - function handle of biased (dimensionless) potential for k-th umbrella-window
%               [cell nwindow x 1]
% * data_kn   - coordinates relevant to biased potentials
%               [cell nwindow x 1]
% * f_k       - (dimensionless) free energies of umbrella-windows
%               [double K x 1]
% 
%% Example
%#
% 
%% See also
% mbarpmf mbarexpectation
%
%% References
% [1] M. R. Shirts and J. D. Chodera, J Chem Phys 129, 124105 (2008).
% [2] Z. Tan, Journal of the American Statistical Association 99, 1027 (2004).
% [3] C. H. Bennett, J Comput Phys 22, 245 (1976).
%

% The names of variables and indicies follow the convention of Ref 1.
% k, l, i, j run over the umbrella-windows.
% n runs over the snapshots.

%% preparation
% K: number of umbrella windows
K = numel(u_kn); 
assert(K == numel(fhandle_k), 'the numbers of umbrella windows in u_kn and fhandle do not match...');
assert(K == numel(data_kn), 'the numbers of umbrella windows in u_kn and data_kn do not match...');

% N_k: number of data in k-th umbrella window
N_k = zeros(K, 1);
for k = 1:K
  N_k(k) = size(u_kn{k}, 1);
end
N_max = max(N_k);

% tolerance to check convergence of iterations
if (nargin < 4) || numel(tolerance) == 0
  tolerance = 10^(-8);
end

%% calculate energies evaluated at different umbrella windows
u_kln = zeros(K, K, N_max);
for k = 1:K
  for l = 1:K
    for n = 1:N_k(k);
      u_kln(k, l, n) = u_kn{k}(n) + fhandle_k{l}(data_kn{k}(n, :));
    end
  end
end

%% solve the MBAR equation by self-consistent iteration
f_k = zeros(K, 1);
f_k_new = f_k;

log_wi_jn = zeros(K, N_max);
for j = 1:K
  log_wi_jn(j, 1:N_k(j)) = 1;
end
index = log_wi_jn > 0;

for count_iteration = 1:5
  for i = 1:K
    log_wi_jn = calc_log_wi_jn(N_k, f_k, u_kln, squeeze(u_kln(:, i, :)), K, N_max);
    f_k_new(i) = - logsumexp(log_wi_jn(index));
  end
  
  f_k_new = f_k_new - f_k_new(1);
  check_convergence = max(abs(f_k_new - f_k))./std(f_k_new);
  f_k = f_k_new;

  fprintf('%dth iteration  delta = %e  tolerance = %e\n', count_iteration, check_convergence, tolerance);
  fprintf('free energies = ');
  fprintf('%f ', f_k);
  fprintf('\n');
  fprintf('\n');
end

%% solve the MBAR equation by the Newton-Raphson method
first_gamma = 0.1;
gamma = 1.0;

%f_k = zeros(K, 1);
f_k_new = f_k;
check_convergence = inf;

log_wi_jn = zeros(K, N_max);
for j = 1:K
  log_wi_jn(j, 1:N_k(j)) = 1;
end
index = log_wi_jn > 0;
W_nk = zeros(sum(N_k), K);

count_iteration = 5;
while check_convergence > tolerance
  f_k_new = f_k;
  for i = 1:K
    log_wi_jn = calc_log_wi_jn(N_k, f_k, u_kln, squeeze(u_kln(:, i, :)), K, N_max);
    W_nk(:,i) = exp(log_wi_jn(index) + f_k(i));
  end

  g = zeros(K-1, 1);   % gradient
  H = zeros(K-1, K-1); % hessian
  for i = 1:(K-1)
    g(i) = N_k(i+1) - N_k(i+1) * sum(W_nk(:,i+1));
    H(i,i) = - sum(N_k(i+1) * W_nk(:,i+1) .* (1.0 - N_k(i+1) * W_nk(:,i+1)));
    for j = 1:(i-1)
      H(i,j) = sum((N_k(i+1) * W_nk(:,i+1)) .* (N_k(j+1) * W_nk(:,j+1)));
      H(j,i) = H(i,j);
    end
  end

  Hinvg = pinv(H, 1.0e-10) * g;
  for k = 1:(K-1)
    if (count_iteration == 0)
      f_k_new(k+1) = f_k_new(k+1) - first_gamma * Hinvg(k);
    else
      f_k_new(k+1) = f_k_new(k+1) - gamma * Hinvg(k);
    end
  end

  check_convergence = max(abs(f_k_new - f_k))./std(f_k_new);
  f_k = f_k_new;

  count_iteration = count_iteration + 1;
  fprintf('%dth iteration  delta = %e  tolerance = %e\n', count_iteration, check_convergence, tolerance);
  fprintf('free energies = ');
  fprintf('%f ', f_k);
  fprintf('\n');
  fprintf('\n');
end

% recompute all free energies
for i = 1:K
  log_wi_jn = calc_log_wi_jn(N_k, f_k, u_kln, squeeze(u_kln(:, i, :)), K, N_max);  
  f_k(i) = - logsumexp(log_wi_jn(index));
end
f_k = f_k - f_k(1);
fprintf('free energies = ');
fprintf('%f ', f_k);
fprintf('\n');
fprintf('\n');


%% calc log weights
function log_wi_jn = calc_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
log_wi_jn = zeros(K, N_max);
for k = 1:K
  log_wi_jn(k, 1:N_k(k)) = - logsumexp2(repmat(log(N_k), 1, N_k(k)) + repmat(f_k, 1, N_k(k)) - (squeeze(u_kln(k, :, 1:N_k(k))) - repmat(u_kn(k, 1:N_k(k)), K, 1)));
end


%% logsumexp (input should be a vector)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;

    
%% logsumexp2 (input should be an array. sums over rows)
function s = logsumexp2(x)
max_x = max(x);
exp_x= exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

