function [f_k, prob_m, center_m, h_km, bias_km, N_m] = wham(edge_m, fhandle_k, data_k, kbt, tolerance)
%% wham
% calculate dimensionless free energies of umbrella windows and the unbiased probabilities on grids by using the weighted histogram analysis method (WHAM)
%
%% Syntax
%# [f_k, prob_m, center_m, h_km, bias_km, N_m] = wham(edge_m, fhandle_k, data_k, kbt)
%# [f_k, prob_m, center_m, h_km, bias_km, N_m] = wham(edge_m, fhandle_k, data_k, kbt, tolerance)
%
%% Description
%
% * edge_m    - edges of bins
%               [double M]
% * fhandle_k - cell of function handles which represent biased potentials
%               [cell K]
% * data_m    - cell of trajectories in a space where histograms are counted
%               [cell K]
% * kbt       - Kb*T in kcal/mol
%               [double scalar]
% * tolerance - tolerance of the convergence in WHAM iterations
%               [double scalar]
% * f_k       - dimensionless free energies of umbrella windows
%               [double K x 1]
% * prob_m    - unbiased probability of bins (prob(m) = exp(-kbt*U(x_m))/Z(kbt) * dx_m)
%               [double 1 x M]
% 
%% Example
%#
% 
%% See also
% mbar
%
%% References
% [1] S. Kumar, D. Bouzida, R. H. Swendsen, P. A. Kollman, and
%     J. M. Rosenberg, J. Comput. Chem. 13, 1011 (1992). 
% [2] B. Roux, Computer Physics Communications 91, 275 (1995).
% [3] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and  
%     K. A. Dill, J. Chem. Theory Comput. 3, 26 (2007).
% [4] J. S. Hub, B. L. de Groot, and D. van der Spoel,
%     J. Chem. Theory Comput. 6, 3713 (2010). 
%

% The names of variables and indicies follow those of Ref [3]. 
% We assume an array structure whose rows correspond to umbrella
% windows and columns are bins. 

%% preparation
% K: number of umbrella windows
K = numel(data_k); 
assert(K == numel(fhandle_k), 'the numbers of windows in data and fhandle do not match...');
% M: number of bins
M = numel(edge_m) - 1;
% calculate the centers of bins
center_m = edge_m + 0.5*(edge_m(2) - edge_m(1));
center_m(end) = [];
% tolerance to check convergence of iterations
if (nargin < 5) | numel(tolerance) == 0
  tolerance = 10^(-8);
end

%% count histogram
h_km = zeros(K, M);
for k = 1:K
  hh = histc(data_k{k}, edge_m);
  h_km(k, :) = hh(1:(end-1))';
end
% check the number of histogram counts for each window
N_m = sum(h_km, 2);

%% calculate bias-energy
bias_km = zeros(K, M);
for k = 1:K
  for m = 1:M
    bias_km(k, m) = fhandle_k{k}(center_m(m));
  end
end
% convert to dimensionless energies
bias_km = bias_km ./ kbt;

%% calculate statistical inefficiency


%% solve the WHAM equations by self-consistent iteration
f_k = zeros(K, 1);
prob_m = zeros(1, M);
check_convergence = inf;

count_iteration = 0;
while check_convergence > tolerance
  eq1_inv_km  = exp(bsxfun(@minus, f_k, bias_km));
  eq1_inv_m   = sum(bsxfun(@times, N_m, eq1_inv_km));
  prob_m  = sum(h_km)./eq1_inv_m;
  f_k_new = - logsumexp(-bias_km' + repmat(log(prob_m'), 1, K));
  f_k_new = f_k_new';
  f_k_new = f_k_new - f_k_new(1);
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

% normalize unbiased probability
prob_m = prob_m./sum(prob_m);

%% bootstrap?
f_k_error = 0;
prob_m_error = 0;


%% logsumexp (input should be array. sums over rows)
function s = logsumexp(x)
[K, M] = size(x);
max_x = max(x);
exp_x= exp(x - repmat(max_x, K, 1));
s = log(sum(exp_x)) + max_x;

