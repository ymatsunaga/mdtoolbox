function pmf_i = mbarpmf(u_kl, bin_k, f_k, u_k)
%% mbarpmf
% calculate the potential of mean force of bins by using the Multistate Bennet Acceptance Ratio Method (MBAR)
%
%% Syntax
%# pmf = mbarpmf(u_kl, bin_k, f_k)
%# pmf = mbarpmf(u_kl, bin_k, f_k, u_k)
%
%% Description
%
% * u_kl      - reduced potential energy of umbrella simulation k evaluated at umbrella l
%               [cell numbrella x numbrella]
% * bin_k     - binned trajectories of umbrella simulation k
%               [cell numbrella x 1]
% * f_k       - reduced free energies of umbrella k
%               [double numbrella x 1]
% * u_k       - reduced potential energy umbrella simulation k under the condition where PMF is evaluated 
%               [cell numbrella x 1]
%               If omitted, u_k = 0 is assumed
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

% The names of variables and indicies follow the convention of Ref 1.
% We assume array structures whose rows correspond to umbrella
% windows and columns are bins. 

%% preparation
% K: number of umbrella windows
[K, L] = size(u_kl); 
assert(K == L, 'the numbers of rows and columns of u_kl should be same...');

% N_k: number of data in k-th umbrella window
N_k = zeros(K, 1);
for k = 1:K
  N_k(k) = size(u_kl{k, 1}, 1);
end
N_max = max(N_k);

% conversion from cell (u_kl) to array (u_kln)
u_kln = zeros(K, K, N_max);
for k = 1:K
  for l = 1:K
    u_kln(k, l, 1:N_k(k)) = u_kl{k, l};
  end
end
clear u_kl;

% conversion from cell (u_k) to array (u_kn)
u_kn = zeros(K, N_max);
if exist('u_k', 'var') && ~isempty(u_k);
  for k = 1:K
    u_kn(k, 1:N_k(k)) = u_k{k};
  end
  clear u_k;
end

% conversion from cell (bin_k) to array (bin_kn)
bin_kn = zeros(K, N_max);
for k = 1:K
  bin_kn(k, 1:N_k(k)) = bin_k{k};
end
clear bin_k;

%% calc PMF
log_w_kn = zeros(K, N_max);
for j = 1:K
  log_w_kn(j, 1:N_k(j)) = 1;
end
index = log_w_kn > 0;

log_w_kn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max);
log_w_n  = log_w_kn(index);

nbin   = max(bin_kn(:));
pmf_i  = zeros(nbin, 1);
for ibin = 1:nbin
  lindex_ibin = (bin_kn(index) == ibin);
  if any(lindex_ibin)
    pmf_i(ibin) = - logsumexp(log_w_n(lindex_ibin));
  else
    pmf_i(ibin) = NaN;
  end
end
%pmf_i = pmf_i - pmf_i(1);


%% logsumexp (input should be vector)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;

