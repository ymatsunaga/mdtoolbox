function [f_k, log_prob_m, center_m, h_km, bias_km, N_k] = wham(edge_m, fhandle_k, data_k, kbt, tolerance)
%% wham
% calculate the dimensionless free energies of umbrella windows and the unbiased probabilities on grids by using the WHAM
%
%% Syntax
%# [f_k, prob_m, center_m, h_km, bias_km, N_k] = wham(edge_m, fhandle_k, data_k, kbt)
%# [f_k, prob_m, center_m, h_km, bias_km, N_k] = wham(edge_m, fhandle_k, data_k, kbt, tolerance)
%
%% Description
%
% * edge_m    - edges of bins
%               [double M]
% * fhandle_k - cell of function handles which represent biased potentials
%               [cell K]
% * data_k    - cell of trajectories in a space where histograms are counted
%               [cell K]
% * kbt       - Kb*T in kcal/mol
%               [double scalar]
% * tolerance - tolerance of the convergence in WHAM iterations
%               [double scalar]
% * f_k       - dimensionless free energies of umbrella windows
%               [double K x 1]
% * log_prob_m  - log of unbiased probability in bins (prob(m) = exp(-kbt*U(x_m))/Z(kbt) * dx_m)
%                 [double 1 x M]
% 
%% Example
%#
% 
%% See also
% ptwham
%
%% References
% [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and  
%     K. A. Dill, J. Chem. Theory Comput. 3, 26 (2007).
% [2] S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and
%     J. M. Rosenberg, J. Comput. Chem. 13, 1011 (1992). 
% [3] B. Roux, Computer Physics Communications 91, 275 (1995).
% [4] J. S. Hub, B. L. de Groot, and D. van der Spoel,
%     J. Chem. Theory Comput. 6, 3713 (2010). 
%

% The notation of variables and indicies follows Ref [1]. 
% Also, we assume an array structure whose 
% rows(k) correspond to umbrella windows and columns(m) are bins. 

%% preparation
% K: number of umbrella windows
K = numel(data_k); 
% M: number of bins
M = numel(edge_m) - 1;
% center_m: centers of data bins (m)
center_m = edge_m + 0.5*(edge_m(2) - edge_m(1));
center_m(end) = [];
% tolerance to check convergence of iterations
if (nargin < 5) || numel(tolerance) == 0
  tolerance = 10^(-8);
end

% check consistency of the number of umbrella-windows
if K ~= numel(fhandle_k)
  error('# of umbrella windows do not match... data has %d windows. fhandle has %d windows.', numel(data_k), numel(fhandle_k));
end

%% calculate histogram (h_km)
% h_km: number of smaples in data-bin (m) from umbrella-window (k)
h_km = zeros(K, M);
for k = 1:K
  h_m = histc(data_k{k}, edge_m);
  h_km(k, :) = h_m(1:(end-1))';
end

%% calculate counts (N_k)
% N_k: number of samples counted in the histogram from umbrella-window (k)
N_k = sum(h_km, 2);
log_N_k = zeros(K, 1);
log_N_k(N_k > 0) = log(N_k(N_k > 0));

%% calculate statistical inefficiency (g_km)
% g_km: statistical efficiency of indicator function for umbrella-window (k) in data-bin (m)
g_km = ones(K, M);

%% calculate effective histogram (heff_m)
% heff_m: effective number of independent samples in data-bin (m)
heff_m = sum(h_km./g_km);
log_heff_m = zeros(1, M);
log_heff_m(:) = double(log(eps('single')));
log_heff_m(heff_m > 0) = log(heff_m(heff_m > 0));

%% calculate effective counts (Neff_km)
% Neff_km: effective number of independent samples for umbrella window (k) in data-bin (m)
Neff_km = bsxfun(@rdivide, N_k, g_km);
log_Neff_km = zeros(K, M);
log_Neff_km(:) = double(log(eps('single')));
log_Neff_km(Neff_km > 0) = log(Neff_km(Neff_km > 0));

%% calculate bias-factor
% bias_km: bias-factor for umbrella window (k) in data-bin (m)
bias_km = zeros(K, M);
for k = 1:K
  for m = 1:M
    bias_km(k, m) = fhandle_k{k}(center_m(m));
  end
end
% convert to dimensionless energies
bias_km = bias_km ./ kbt;

%% solve the WHAM equations by self-consistent iteration
f_k = zeros(K, 1);
log_prob_m = zeros(1, M);
check_convergence = inf;

count_iteration = 0;
while check_convergence > tolerance

  % 1st equation
  log_denominator_km = bsxfun(@plus, f_k, -bias_km);
  log_denominator_km = log_denominator_km + log_Neff_km;
  log_denominator_m  = logsumexp(log_denominator_km);

  log_numerator_m = log_heff_m;

  log_prob_m = log_numerator_m - log_denominator_m;
  
  % 2nd equation
  f_k_new = - logsumexp(bsxfun(@plus, log_prob_m', -bias_km'));
  f_k_new = f_k_new';
  f_k_new = f_k_new - f_k_new(1);

  % check convergence
  check_convergence = max(abs(f_k_new - f_k))./std(f_k_new);
  f_k = f_k_new;

  count_iteration = count_iteration + 1;
  if mod(count_iteration, 100) == 0
    fprintf('%dth iteration\n', count_iteration);
    fprintf('free energies = ');
    fprintf('%f ', f_k);
    fprintf('\n');
  end
end

%% bootstrap?
% f_k_error = 0;
% prob_m_error = 0;


%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

