function pmf_i = mbarpmf(u_kn, fhandle_k, data_kn, bin_kn, f_k)
%% mbarpmf
% calculate the potential of mean force of bins
%
%% Syntax
%# pmf = mbarpmf(u_kn, fhandle_k, data_kn, bin_kn, f_k)
%
%% Description
%
% * u_kn      - u (x_{kn})
%               [cell nwindow x nstep]
% * fhandle   - cell of function handles which represent biased potentials
%               [cell nwindow]
% * bin       - cell of trajectories in a space where histograms are counted
%               [cell nwindow]
% 
%% Example
%#
% 
%% See also
% mbar mbarexpectation
%
%% References
% [1] M. R. Shirts and J. D. Chodera, J Chem Phys 129, 124105 (2008).
% [2] Z. Tan, Journal of the American Statistical Association 99, 1027 (2004).
% [3] C. H. Bennett, Journal of Computational Physics (1976).
%

% The names of variables and indicies follow those of Ref [1]. 
% We assume array structures whose rows correspond to umbrella
% windows and columns are bins. 

%% preparation
% K: number of umbrella windows
K = numel(u_kn); 
assert(K == numel(fhandle_k), 'the numbers of umbrella windows in u_kn and fhandle do not match...');
assert(K == numel(bin_kn), 'the numbers of umbrella windows in u_kn and fhandle do not match...');

% N_k: number of data in k-th umbrella window
N_k = zeros(K, 1);
for k = 1:K
  N_k(k) = size(u_kn{k}, 1);
end
N_max = max(N_k);

u_cell_kn = u_kn;
u_kn = zeros(K, N_max);
for k = 1:K
  for n = 1:N_k(k);
    u_kn(k, n) = u_cell_kn{k}(n);
  end
end
clear u_cell_kn;

bin_cell_kn = bin_kn;
bin_kn = zeros(K, N_max);
for k = 1:K
  for n = 1:N_k(k);
    bin_kn(k, n) = bin_cell_kn{k}(n);
  end
end
clear bin_cell_kn;

%% calculate energies evaluated at different umbrella windows
u_kln = zeros(K, K, N_max);
for k = 1:K
  for l = 1:K
    for n = 1:N_k(k);
      u_kln(k, l, n) = u_kn(k, n) + fhandle_k{l}(data_kn{k}(n, :));
    end
  end
end

%% calc PMF
log_w_kn = zeros(K, N_max);
for j = 1:K
  log_w_kn(j, 1:N_k(j)) = 1;
end
index = log_w_kn > 0;

log_w_kn = calc_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max);
log_w_n  = log_w_kn(index);

nbin   = max(bin_kn(:));
pmf_i  = zeros(nbin, 1);
for ibin = 1:nbin
  lindex_ibin = (bin_kn(index) == ibin);
  pmf_i(ibin) = - logsumexp(log_w_n(lindex_ibin));
end
pmf_i = pmf_i - pmf_i(1);


%% calc log weights
function log_wi_jn = calc_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
log_wi_jn = zeros(K, N_max);
for k = 1:K
  log_wi_jn(k, 1:N_k(k)) = - logsumexp2(repmat(log(N_k), 1, N_k(k)) + repmat(f_k, 1, N_k(k)) - (squeeze(u_kln(k, :, 1:N_k(k))) - repmat(u_kn(k, 1:N_k(k)), K, 1)));
end


%% logsumexp (input should be vector)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;

    
%% logsumexp2 (input should be array. sum over rows)
function s = logsumexp2(x)
[K, M] = size(x);
max_x = max(x);
exp_x = exp(x - repmat(max_x, K, 1));
s = log(sum(exp_x)) + max_x;

