function [f_k, prob_m, center_m, h_km, bias_km, N_m] = calcwham(data_k, fhandle_k, edge_m, kbt, tolerance)
%% calcwham
% calculate the free energy differences of umbrella windows and the unbiased probabilities on grids using the weighted histogram analysis method (WHAM)
%
%% Syntax
%# [f, prob, center] = calcwham(data, fhandle, kbt, edge)
%
%% Description
%
% * data      - unbiased potential energy [cell nwindow]
%               [cell nwindow]
% * fhandle   - cell of function handles of biased potentials
%               [cell nwindow]
% * edge      - edges of bins
%               [double nbin]
% * kbt       - Kb*T in kcal/mol
%               [double scalar]
% * tolerance - tolerance of the convergence in WHAM iterations
%               [double scalar]
% * f         - free energy constants of umbrella windows
%               [double nwindow x 1]
% * prob      - unbiased probability of bins
%               [double nbin x 1]
% 
%% Example
%#
% 
%% See also
% calcmbar
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
% We assume array structures whose rows correspond to umbrella windows and columns are bins

%% preparation
% number of windows
K = numel(data_k); 
assert(K == numel(fhandle_k), 'the numbers of windows in data and fhandle do not match...');
% number of bins
M = numel(edge_m) - 1;
% calc the centers of bins
center_m = edge_m + 0.5*(edge_m(2) - edge_m(1));
center_m(end) = [];
% temperature
beta = 1./kbt;
% tolerance to check convergence of iterations
if (nargin < 5) | numel(tolerance) == 0
  tolerance = 10^(-10);
end

%% calculate histogram
h_km = zeros(K, M);
for k = 1:K
  hh = histc(data_k{k}, edge_m);
  h_km(k, :) = hh(1:(end-1))';
end

%% calculate the number of data used in the histogram for each window
N_m = sum(h_km, 2);

%% calculate bias-energy
bias_km = zeros(K, M);
for k = 1:K
  for m = 1:M
    bias_km(k, m) = fhandle_k{k}(center_m(m));
  end
end

%% solve the WHAM equations by iterations
f_k = zeros(K, 1);
prob_m = zeros(1, M);
check_convergence = inf;

while check_convergence > tolerance
  x_inv   = exp(beta*bsxfun(@minus, f_k, bias_km));
  x_inv   = sum(bsxfun(@times, N_m, x_inv));
  prob_m  = sum(h_km)./x_inv;
  f_k_new = - kbt*logsumexp2(-beta*bias_km' + repmat(log(prob_m'), 1, K));
  %f_k_new = sum(exp(-beta*bias_km').*repmat(log(prob_m'), 1, K));
  %f_k_new = - kbt*log(f_k_new);
  f_k_new = f_k_new - f_k(1)
  f_k_new = f_k_new';
  check_convergence = max(abs(f_k_new - f_k))./std(f_k_new);
  f_k = f_k_new;
end

% normalize unbiased probability
prob_m = prob_m./sum(prob_m);

%% bootstrap?
f_k_error = 0;
prob_m_error = 0;


%% logsumexp
function s = logsumexp(log_terms)
max_log_term = max(log_terms);
terms = exp(log_terms - max_log_term);
s = log(sum(terms)) + max_log_term;

    
%% logsumexp2 (a bit faster than logsumexp)
function s = logsumexp2(log_terms)
[K, M] = size(log_terms);
max_log_term = max(log_terms);
terms = exp(log_terms - repmat(max_log_term, K, 1));
s = log(sum(terms)) + max_log_term;


