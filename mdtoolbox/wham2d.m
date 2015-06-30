function [f_k, pmf, center_mx, center_my] = wham2d(data_k, fhandle_k, temperature, edge_mx, edge_my)
%% wham2d
% calculate 2-dimensional (dimensionless relative) free energies and potential of mean force in 2-dimensional data-bins by using the WHAM
%
%% Syntax
%# [f_k, pmf, center_mx, center_my] = wham2d(data_k, fhandle_k, temperature, edge_mx, edge_my)
%
%% Description
%
% * data_k      - cell of trajectories in a space where histograms are counted
%                 [cell K x 1]
% * fhandle_k   - cell of function handles which represent biased potentials in the x-axis
%                 [cell K x 1]
% * temperature - Temperature in Kelvin
%                 [double scalar]
% * edge_mx     - edges of data-bins in the x-axis
%                 [double 1 x MX or MX x 1]
% * edge_my     - edges of data-bins in the y-axis
%                 [double 1 x MY or MY x 1]
% 
% * f_k         - dimensionless free energies of umbrella-windows in 2-dimensional space
%                 [double K x 1]
% * pmf         - potential of mean force in data-bins under unbiased condition in 2-dimensional space
%                 [double MY x MX]
% * center_mx   - centers of data-bins in the x-axis
%                 [double MX x 1]
% * center_my   - centers of data-bins in the y-axis
%                 [double MY x 1]
% 
%% Example
%# see http://mdtoolbox.readthedocs.org/en/latest/alad_2D_umbrella_wham.html
% 
%% See also
% wham ptwham mbar mbarpmf
%
%% References
% This function uses the method described in
% [1] S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and
%     J. M. Rosenberg, J. Comput. Chem. 13, 1011 (1992). 
% Also, some parts are based on the following papers:
% [2] B. Roux, Computer Physics Communications 91, 275 (1995).
% [3] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and  
%     K. A. Dill, J. Chem. Theory Comput. 3, 26 (2007).
% [4] J. S. Hub, B. L. de Groot, and D. van der Spoel,
%     J. Chem. Theory Comput. 6, 3713 (2010). 
%
%% TODO
% bootstrap
%

% The names of variables and indicies follow the convention of Ref 1.
% Also, we assume an array structure whose 
% rows(k) correspond to umbrella-windows and columns(m) are bins. 

%% preparation
% Boltzmann constant in kcal/(mol K)
C = getconstants();
KB = C.KB;
% Tolerance for the convergence of iterations
TOLERANCE = 10^(-8);

% convert from array to cell
if ~iscell(data_k)
  K = size(data_k, 1);
  t = data_k;
  data_k = cell(K, 1);
  for k = 1:size(data_k, 1)
    data_k{k} = t(k, :);
  end
  clear t;
end

% K: number of umbrella-windows
K = numel(data_k);

% check consistency of the number of umbrella-windows
K2 = numel(fhandle_k);
assert(K == K2, sprintf('# of umbrella-windows do not match... data has %d windows. fhandle has %d windows.', K, K2));

% MX: number of data-bins in the x-axis
if ~exist('edge_mx', 'var') || isempty(edge_mx)
  MX = 50;
  edge_min = zeros(K, 1);
  edge_max = zeros(K, 1);
  for k = 1:K
    edge_min(k) = min(data_k{k}(:, 1));
    edge_max(k) = max(data_k{k}(:, 1));
  end
  edge_min = min(edge_min);
  edge_max = max(edge_max);
  edge_mx = linspace(edge_min, edge_max+(MX*eps), MX+1)';
else
  MX = numel(edge_mx) - 1;
end
center_mx = 0.5 * (edge_mx(2:end) + edge_mx(1:(end-1)));

% MY: number of data-bins in the y-axis
if ~exist('edge_my', 'var') || isempty(edge_my)
  MY = 50;
  edge_min = zeros(K, 1);
  edge_max = zeros(K, 1);
  for k = 1:K
    edge_min(k) = min(data_k{k}(:, 2));
    edge_max(k) = max(data_k{k}(:, 2));
  end
  edge_min = min(edge_min);
  edge_max = max(edge_max);
  edge_my = linspace(edge_min, edge_max+(MY*eps), MY+1)';
else
  MY = numel(edge_my) - 1;
end
center_my = 0.5 * (edge_my(2:end) + edge_my(1:(end-1)));

% M: number of data-bins
M = MX * MY;

%% calculate histogram (h_km)
% h_km: number of samples in data-bin (m) from umbrella-window (k)
h_km = zeros(K, M);
for k = 1:K
  h_mx_my = hist3(data_k{k}, {edge_mx, edge_my});
  h = h_mx_my(1:MX, 1:MY);
  h_km(k, :) = h(:)';
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
% Neff_km: effective number of independent samples for umbrella-window (k) in data-bin (m)
Neff_km = bsxfun(@rdivide, N_k, g_km);
log_Neff_km = zeros(K, M);
log_Neff_km(:) = double(log(eps('single')));
log_Neff_km(Neff_km > 0) = log(Neff_km(Neff_km > 0));

%% calculate bias-factor
% bias_km: bias-factor for umbrella-window (k) in data-bin (m)
bias_km = zeros(K, M);
for k = 1:K
  for mx = 1:MX
    for my = 1:MY
      m = (my-1)*MX + mx;
      bias_km(k, m) = fhandle_k{k}([center_mx(mx) center_my(my)]);
    end
  end
end
% convert to dimensionless energies
bias_km = bias_km ./ (KB*temperature);

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

pmf = reshape(log_prob_m, MX, MY);
pmf = - KB * temperature * pmf;
pmf = pmf - min(pmf(:));
pmf = pmf';

%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

