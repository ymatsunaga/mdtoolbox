function [pmf_j, center_j, w_kn] = ptwhampmf(data_kn, temperature, log_state_m, center_m, h_km, bin_kn)
%% ptwhampmf
% calculate potential of mean force using the results of PTWHAM
%
%% Syntax
%# [pmf_j, center_j, w_kn] = ptwhampmf(data_kn, temperature, log_state_m, center_m, h_km, bin_kn)
%
%% Description
%
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
% rows correspond to umbrella-windows and columns are energy-bins. 

%% preparation
% Boltzmann constant in kcal/(mol K)
KB = 0.00198719168260038;
K = numel(data_kn);
M = numel(log_state_m);
temp0 = temperature;
beta0 = 1./(KB*temp0);

%% accumulate histogram
h_m = sum(h_km, 1)';
log_h_m = zeros(M, 1);
log_h_m(:) = double(log(eps('single')));
log_h_m(h_m > 0) = log(h_m(h_m > 0));

%% calculate weights
% w_kn: weight of the n-the snapshot from umbrella-window (k)
for k = 1:K
  N = numel(data_kn{k});
  w_n = zeros(N, 1);
  for n = 1:N
    m = bin_kn{k}(n);
    w_n(n) = log_state_m(m) - log_h_m(m) - beta0*center_m(m);
  end
  w_kn{k} = w_n;
end

% normalize weights by logsumexp
w_kn_max = max(cellfun(@max, w_kn));
for k = 1:K
  w_kn{k} = exp(w_kn{k} - w_kn_max);
end
w_kn_sum = sum(cellfun(@sum, w_kn));
for k = 1:K
  w_kn{k} = w_kn{k} ./ w_kn_sum;
end

%% calculate PMF
J = 100;
edge_min = min(cellfun(@min, data_kn));
edge_max = max(cellfun(@max, data_kn));
edge_j = linspace(edge_min, edge_max+(J*eps), J+1)';
center_j = 0.5 * (edge_j(2:end) + edge_j(1:(end-1)));

for k = 1:K
  [~, bin_n] = histc(data_kn{k}, edge_j);
  bindata_kn{k} = bin_n;
end

pmf_j = zeros(J, 1);
for k = 1:K
  N = numel(data_kn{k});
  for n = 1:N
    j = bindata_kn{k}(n);
    pmf_j(j) = pmf_j(j) + w_kn{k}(n);
  end
end

pmf_j = - KB*temp0*log(pmf_j);


