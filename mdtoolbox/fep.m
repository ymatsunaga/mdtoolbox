function f_k = fep(delta_u, T, nblock)
%% fep
% compute free energy difference by free energy perturbation (exponential averaging)
%
%% Syntax
%# f_k = fep(delta_u, T)
%
%% Description
% free energy difference is estimated by exponential averaging
% error is estimated by standard point estimation theory
% 
%% Example
% 
%% See also
%
%% References
%
%

%%%%%%%% main
if ~exist('nblock', 'var') || isempty(nblock)
  nblock = 1;
end

K = numel(delta_u);

s = getconstants();
KBT = s.KB*T;
beta = 1./KBT;

% EXP for each block
f_k = zeros(K, 2);
for k = 1:K
  if nblock > 1
    [f, f_var] = kernelfunction(delta_u{k}, beta, nblock);
  else
    [f, f_var] = kernelfunction_noblock(delta_u{k}, beta);
  end
  if k ==1
    f_k(k, 1) = f;
    f_k(k, 2) = f_var;
  else
    f_k(k, 1) = f_k(k-1, 1) + f;
    f_k(k, 2) = f_k(k-1, 2) + f_var;
  end
end

f_k(:, 2) = sqrt(f_k(:, 2));

%%%%%%%% kernel function
function [f, f_var] = kernelfunction_noblock(u, beta)
nstep = numel(u);
w = - beta*u;
f = - (1./beta) * logsumexp(w) + (1./beta)*log(nstep);
w_mean = mean(exp(w));
w_var = var(exp(w));
f_var = (1./nstep)*w_var./(w_mean.^2);

function [f, f_var] = kernelfunction(u, beta, nblock)
nstep = numel(u);
interface = round(linspace(0, nstep, nblock+1));
index_n = {};
for n = 1:nblock
  index_n{n} = (interface(n)+1):interface(n+1);
end

f_n = [];
for n = 1:nblock
  nstep = numel(index_n{n});
  w = - beta*u(index_n{n});
  f_n = [f_n; - (1./beta) * logsumexp(w) + (1./beta)*log(nstep)];
end
f = mean(f_n);
f_var = var(f_n);

%%%%%%%% logsumexp (input should be a vector)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;

