function log_wi_jn = mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, K, N_max)
%% mbar_log_wi_jn
% An auxiliary function for mbar routines (mbar*.m)
%
%% References
% [1] M. R. Shirts and J. D. Chodera, J Chem Phys 129, 124105 (2008).
% [2] Z. Tan, Journal of the American Statistical Association 99, 1027 (2004).
% [3] C. H. Bennett, Journal of Computational Physics (1976).
%

% The names of variables and indicies follow the convention of Ref 1.
% We assume array structures whose rows correspond to umbrella
% windows and columns are bins. 

log_wi_jn = zeros(K, N_max);
for k = 1:K
  log_wi_jn(k, 1:N_k(k)) = - logsumexp2(repmat(log(N_k), 1, N_k(k)) + repmat(f_k, 1, N_k(k)) - (squeeze(u_kln(k, :, 1:N_k(k))) - repmat(u_kn(k, 1:N_k(k)), K, 1)));
end

%% logsumexp2 (input should be an array. sums over rows)
function s = logsumexp2(x)
max_x = max(x);
exp_x= exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

