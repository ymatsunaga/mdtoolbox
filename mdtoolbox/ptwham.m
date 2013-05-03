function [f_l, log_state_m, center_m, h_km, bias_lm, N_kl] = ptwham(edge_m, temperature_k, penergy_k)
%% ptwham
% calculate dimensionless free energies of umbrella-windows and the potential energy density of states by using the PTWHAM
%
%% Syntax
%# [f_l, log_state_m, center_m, h_km, bias_lm, N_kl] = wham(edge_m, temperature_k, penergy_k)
%# [f_l, log_state_m, center_m, h_km, bias_lm, N_kl] = wham(edge_m, temperature_k, penergy_k, tolerance)
%
%% Description
%
% * edge_m        - edges of energy-bins
%                   [double M]
% * temperature_k - cell of trajectories in a space where histograms are counted
%                   [cell K]
% * data_k        - cell of trajectories in a space where histograms are counted
%                   [cell K]
% * kbt           - Kb*T in kcal/mol
%                   [double scalar]
% * f_l           - dimensionless free energies of umbrella-windows
%                   [double K x 1]
% * log_state_m    - log of unbiased probability in energy-bins (prob(m) = exp(-kbt*U(x_m))/Z(kbt) * dx_m)
%                    [double 1 x M]
% 
%% Example
%#
% 
%% See also
% wham
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
% rows correspond to umbrella-windows and columns are energy-bins. 

%% preparation
% Boltzmann constant in kcal/(mol K)
KB = 0.00198719168260038;
% K: number of umbrella-windows
K = numel(temperature_k); 
% M: number of energy-bins
M = numel(edge_m) - 1;
% calculate the centers of energy-bins
edge_diff = edge_m(2) - edge_m(1);
center_m = edge_m + 0.5*edge_diff;
center_m(end) = [];
% tolerance for the convergence of iterations
tolerance = 10^(-8);

% check consistency of data sizes
if K ~= numel(penergy_k)
  error('# of umbrella-windows do not match... temperature has %d windows. penergy has %d windows.', numel(data_k), numel(fhandle_k));
end

%% extract unique temperatures
% L: number of unique target temperatures (temp0_l)
temp0_l = [];
for k = 1:K
  temp0_l = unique([temp0_l; temperature_k{k}]);
end
temp0_l = sort(unique(temp0_l(:)))
L = numel(temp0_l);
beta0_l = 1.0./(temp0_l*KB);

%% calculate histogram (h_km)
% h_km: number of smaples in energy-bin (m) from umbrella-window (k)
h_km = zeros(K, M);
for k = 1:K
  h_m = histc(penergy_k{k}, edge_m);
  h_km(k, :) = h_m(1:(end-1))';
end

%% calculate counts (N_kl)
% N_kl: number of temperature temp0_l visited by umbrella-window (k)
N_kl = zeros(K, L);
for k = 1:K
  for l = 1:L
    N_kl(k, l) = sum(temperature_k{k} == temp0_l(l));
  end
end

%% calculate statistical inefficiency (g_km)
% g_km: statistical efficiency of indicator function for umbrella-window (k) in data-bin (m)
g_km = ones(K, M);

%% calculate effective histogram (heff_m)
% heff_m: effective number of independent samples in data-bin (m)
heff_m = sum(h_km./g_km);
log_heff_m = zeros(1, M);
log_heff_m(:) = double(log(eps('single')));
log_heff_m(heff_m > 0) = log(heff_m(heff_m > 0));

%% calculate effective counts (Neff_lm)
% Neff_lm: effective number of independent samples in data-bin (m) at temp0_l
Neff_lm = zeros(L, M);
for l = 1:L
  for m = 1:M
    Neff_lm(l, m) = sum(N_kl(:, l) .* g_km(:, m));
  end
end
log_Neff_lm = zeros(L, M);
log_Neff_lm(:) = double(log(eps('single')));
log_Neff_lm(Neff_lm > 0) = log(Neff_lm(Neff_lm > 0));

%% calculate bias-factor
bias_lm = zeros(L, M);
for l = 1:L
  for m = 1:M
    bias_lm(l, m) = beta0_l(l)*center_m(m);
  end
end

%% solve the WHAM equations by self-consistent iteration
f_l = ones(L, 1);
check_convergence = inf;

count_iteration = 0;
while check_convergence > tolerance

  % 1st equation
  log_denominator_lm = bsxfun(@plus, f_l, -bias_lm);
  log_denominator_lm = log_denominator_lm + log_Neff_lm;
  log_denominator_m  = logsumexp(log_denominator_lm);

  log_numerator_m = log_heff_m;

  log_state_m = log_numerator_m - log_denominator_m;
  
  % 2nd equation
  f_l_new = - logsumexp(bsxfun(@plus, log_state_m', -bias_lm'));
  f_l_new = f_l_new';
  f_l_new = f_l_new - f_l_new(1);

  % check convergence
  check_convergence = max(abs(f_l_new - f_l))./std(f_l_new);
  f_l = f_l_new;

  count_iteration = count_iteration + 1;
  if mod(count_iteration, 1) == 0
    fprintf('%dth iteration\n', count_iteration);
    fprintf('free energies = ');
    fprintf('%f ', f_l);
    fprintf('\n');
  end
end

%% bootstrap?
% f_l_error = 0;
% state_m_error = 0;


%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

