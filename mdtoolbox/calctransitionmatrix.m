function [t, eq, x] = calctransitionmatrix(indexOfCluster, tau, prior)
%% calctransitionmatrix
% calculate transition probability matrix
%
%% Syntax
%# t = calctransitionmatrix(indexOfCluster, tau)
%# t = calctransitionmatrix(indexOfCluster, tau, prior)
%# [t, eq, x] = calctransitionmatrix(indexOfCluster, tau, prior)
%
%% Description
% calculate the transition probability matrix T_ij that transition
% occurs from state i to state j during time tau. 
%
%% Example
%#
% 
%% See also
%
%% References
% [1] K. A. Beauchamp, G. R. Bowman, T. J. Lane, 
%     L. Maibaum, I. S. Haque, and V. S. Pande, 
%     J. Chem. Theory Comput. 7, 3412 (2011).
% 

%% setup
nstep  = numel(indexOfCluster);
nstate = max(indexOfCluster);

index_from = 1:(nstep-tau);
index_to   = (1+tau):nstep;
indexOfCluster_from = indexOfCluster(index_from);
indexOfCluster_to   = indexOfCluster(index_to);

if (nargin <= 2) | (numel(prior) == 0)
  prior = 0;
end

%% calc count matrix
% count transitions and make count matrix C_ij by using a sparse matrix
c = sparse(indexOfCluster_from, indexOfCluster_to, 1, nstate, nstate);

% some regularization (see the references)
s = (c + c') > 0;
c = c + prior*s;
%c = 0.5*(c + c') + prior*s;

%% estimate the symmetric count matrix X_ij using the reversible Maximum Likelihood Estimator
c_sum     = sum(c, 2);
x_sum     = c_sum;
x_sum_old = c_sum;
c_diag    = diag(c);
x_diag    = c_diag;
x         = c;

Q_sym = (c + c');
while true
  rdivide_2ndterm = (nstate*prior + c_sum')./x_sum';
  for i = 1:nstate
    rdivide = (nstate*prior + c_sum(i))./x_sum(i) + rdivide_2ndterm;
    Q = Q_sym(i, :)./rdivide;
    R = 2*prior./rdivide;
    x(i, :) = Q + R;
  end
  x_sum = sum(x, 2);
  x_diag = diag(x);
  if max(abs(x_sum_old - x_sum)) < 10^(-10)
    break
  else
    x_sum_old = x_sum;
  end
end

%% normalize the count matrix [x]
x = x ./ sum(x(:));

%% calc equilibrium probability [eq]
eq = full(sum(x, 2));

%% estimate transition probability matrix [t]
t = x;
for i = 1:nstate
  t(i, :) = x(i, :)./eq(i);
end

