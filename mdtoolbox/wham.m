function [f_k, pmf_m] = wham(h_km, bias_km, g_km)
%% wham
% calculate (reduced) free energies of umbrella-windows and potential of mean force in data-bins by using the WHAM
%
%% Syntax
%# [f_k, pmf_m] = wham(histogram_km, bias_km)
%
%% Description
%
% * histogram_km       - histogram (data counts) of k-th umbrella data counts in m-th bin
%                        [K x M array where K is # of umbrellas and M is # of bins]
% * bias_km (optional) - bias-factor of k-th umbrella-window evaluated at m-th bin-center (DEFAULT = zeros(K, M))
%                        [K x M array where K is # of umbrellas and M is # of bins]
% * g_km (optional)    - statistical efficiency of indicator function for k-th umbrella data in m-th bin (DEFAULT = ones(K, M))
%                        [double 1 x M or M x 1]
%
% * f_k      - reduced relative free energies of umbrella-windows
%              [double K x 1]
% * pmf_m    - reduced potential of mean force in data-bins under unbiased condition
%              [double 1 x M]
% 
%% Example
%# see http://mdtoolbox.readthedocs.org/en/latest/alat_1D_umbrella_wham.html
% 
%% See also
% assign1dbin wham2d ptwham mbar mbarpmf
%
%% References
% This function uses the method described in
% [1] S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and
%     J. M. Rosenberg, J. Comput. Chem. 13, 1011 (1992). 
%

% The names of variables and indicies follow the convention of Ref 1.
% Also, we assume an array structure whose 
% rows(k) correspond to umbrella-windows and columns(m) are bins. 

%% Tolerance for the convergence of iterations
TOLERANCE = 10^(-8);

%% K: number of umbrella-windows
K = size(h_km, 1);

%% M: number of bins
M = size(h_km, 2);

%% check bias_km
if exist('bias_km', 'var') && ~isempty(bias_km)
  assert(size(bias_km, 1) == K, 'sizes of h_km and bias_km are different');
  assert(size(bias_km, 2) == M, 'sizes of h_km and bias_km are different');
else
  bias_km = zeros(K, M);
end

%% check g_km
if exist('g_km', 'var') && ~isempty(g_km)
  assert(size(g_km, 1) == K, 'sizes of h_km and g_km are different');
  assert(size(g_km, 2) == M, 'sizes of h_km and g_km are different');
else
  g_km = ones(K, M);
end

%% calculate counts (N_k)
% N_k: number of samples counted in the histogram from umbrella-window (k)
N_k = sum(h_km, 2);
log_N_k = zeros(K, 1);
log_N_k(N_k > 0) = log(N_k(N_k > 0));

%% calculate effective histogram (heff_m)
% heff_m: effective number of independent samples in data-bin (m)
heff_m = sum(h_km./g_km);
log_heff_m = zeros(1, M);
log_heff_m(:) = double(log(eps('single')));
log_heff_m(heff_m > 0) = log(heff_m(heff_m > 0));

%% calculate effective counts (Neff_km)
% Neff_km: effective number of independent samples for umbrella-window (k) in data-bin (m)
Neff_km = bsxfun(@rdivide, N_k, g_km);
log_Neff_km = zeros(K, M);
log_Neff_km(:) = double(log(eps('single')));
log_Neff_km(Neff_km > 0) = log(Neff_km(Neff_km > 0));

%% solve the WHAM equations by self-consistent iteration
f_k = zeros(K, 1);
check_convergence = inf;

count_iteration = 0;
while check_convergence > TOLERANCE

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
    fprintf('%dth iteration  delta = %e  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
    fprintf('free energies = ');
    fprintf('%f ', f_k);
    fprintf('\n');
    fprintf('\n');
  end
end

pmf_m = - log_prob_m;


%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

