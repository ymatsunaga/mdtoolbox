function [f_k, pmf_m, center_m, h_km, bias_km] = wham(data_k, fhandle_k, temperature, edge_m, nblock)
%% wham
% calculate (dimensionless relative) free energies of umbrella-windows and potential of mean force in data-bins by using the WHAM
%
%% Syntax
%# [f_k, pmf_m, center_m] = wham(data_k, fhandle_k, temperature)
%# [f_k, pmf_m, center_m, h_km, bias_km] = wham(data_k, fhandle_k, temperature, edge_m)
%
%% Description
%
% * data_k      - cell of trajectories in a space where histograms are counted
%                 [cell K x 1]
% * fhandle_k   - cell of function handles which represent biased potentials
%                 [cell K x 1]
% * temperature - Temperature in Kelvin
%                 [double scalar]
% * edge_m      - edges of data-bins
%                 [double 1 x M or M x 1]
% * nblock      - the number of blocks used for error evaluation. default is 1 (no error estimation)
%
% * f_k         - dimensionless relative free energies of umbrella-windows
%                 [double K x 1]
% * pmf_m       - potential of mean force in data-bins under unbiased condition
%                 [double 1 x M]
% 
%% Example
%# see http://mdtoolbox.readthedocs.org/en/latest/alat_1D_umbrella_wham.html
% 
%% See also
% wham2d ptwham mbar mbarpmf
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

%%%%%%%% main
if ~exist('nblock', 'var') || isempty(nblock)
  nblock = 1;
end

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

% M: number of data-bins
if ~exist('edge_m', 'var') || isempty(edge_m)
  M = 100;
  edge_min = min(cellfun(@min, data_k));
  edge_max = max(cellfun(@max, data_k));
  edge_m = linspace(edge_min, edge_max+(M*eps), M+1)';
else
  M = numel(edge_m) - 1;
end
center_m = 0.5 * (edge_m(2:end) + edge_m(1:(end-1)));

% N_k: number of samplec in k-th umbrella window
N_k = zeros(K, 1);
for k = 1:K
  N_k(k) = numel(data_k{k});
end

% WHAM for each block
index_k = {};
f_k = {};
pmf_m = {};
h_km = {};
bias_km = {};
for iblock = 1:nblock
  for k = 1:K
    interface = round(linspace(0, N_k(k), nblock+1));
    index_k{k} = (interface(iblock)+1):interface(iblock+1);
  end
  if nblock > 1
    [f_k{iblock}, pmf_m{iblock}, h_km{iblock}, bias_km{iblock}] ...
        = kernelfunction(data_k, fhandle_k, temperature, edge_m, center_m, K, M, index_k, iblock);
  else
    [f_k{iblock}, pmf_m{iblock}, h_km{iblock}, bias_km{iblock}] ...
        = kernelfunction(data_k, fhandle_k, temperature, edge_m, center_m, K, M, index_k);
  end
  pmf_m{iblock} = pmf_m{iblock}';
end
f_k = cell2mat(f_k);
pmf_m = cell2mat(pmf_m);

if nargout > 1
  if nblock > 1
    [~, index_min] = min(pmf_m(:, 1));
    pmf_m = bsxfun(@minus, pmf_m, pmf_m(index_min, :));
    pmf_m = [mean(pmf_m, 2) 2*std(pmf_m, [], 2)];
  else
    pmf_m = pmf_m - min(pmf_m);
  end
end


%%%%%%%% kernel function
function [f_k, pmf_m, h_km, bias_km] = kernelfunction(data_k, fhandle_k, temperature, edge_m, center_m, K, M, index_k, iblock);

%% Boltzmann constant in kcal/(mol K)
C = getconstants();
KB = C.KB;

%% Tolerance for the convergence of iterations
TOLERANCE = 10^(-8);

%% calculate histogram (h_km)
% h_km: number of samples in data-bin (m) from umbrella-window (k)
h_km = zeros(K, M);
for k = 1:K
  h_m = histc(data_k{k}(index_k{k}), edge_m);
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
% Neff_km: effective number of independent samples for umbrella-window (k) in data-bin (m)
Neff_km = bsxfun(@rdivide, N_k, g_km);
log_Neff_km = zeros(K, M);
log_Neff_km(:) = double(log(eps('single')));
log_Neff_km(Neff_km > 0) = log(Neff_km(Neff_km > 0));

%% calculate bias-factor
% bias_km: bias-factor for umbrella-window (k) in data-bin (m)
bias_km = zeros(K, M);
for k = 1:K
  for m = 1:M
    bias_km(k, m) = fhandle_k{k}(center_m(m));
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
    if exist('iblock', 'var')
      fprintf('[block %d] %dth iteration  delta = %e  tolerance = %e\n', iblock, count_iteration, check_convergence, TOLERANCE);
    else
      fprintf('%dth iteration  delta = %e  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
    end
    fprintf('free energies = ');
    fprintf('%f ', f_k);
    fprintf('\n');
    fprintf('\n');
  end
end

pmf_m = - KB * temperature * log_prob_m;


%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

