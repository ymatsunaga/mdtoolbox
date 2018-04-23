function [t, pi_i] = msmtransitionmatrix(c, TOLERANCE)
%% msmtransitionmatrix
% estimate transition probability matrix from count matrix
%
%% Syntax
%# t = msmtransitionmatrix(c);
%
%% Description
% this routines uses the reversible maximum likelihood estimator
%
%% Example
%#
% 
%% See also
%
%% TODO
% sparse matrix
% prior
%

%% setup
if issparse(c)
  c = full(c);
end

nstate = size(c, 1);

c_sym  = c + c';
x      = c_sym;

c_rs   = sum(c, 2);
x_rs   = sum(x, 2);

if ~exist('TOLERANCE', 'var')
  TOLERANCE = 10^(-4);
end

%% optimization by self-consistent iteration
logL_old = 2*TOLERANCE;
logL = 0;
count_iteration = 0;
x_new = zeros(nstate, nstate);

while abs(logL_old - logL) >= TOLERANCE
  count_iteration = count_iteration + 1;
  logL_old = logL;

  % fixed-point method
  for i = 1:nstate
    for j = 1:nstate
      denom = (c_rs(i)./x_rs(i)) + (c_rs(j)./x_rs(j));
      x_new(i, j) = (c(i, j) + c(j, i)) ./ denom;
    end
  end
  
  % update
  x_rs = sum(x_new, 2);
  x = x_new;
  logL = 0;
  for i = 1:nstate
    for j = 1:nstate
      logL = logL + c(i,j) * log(x(i,j) / x_rs(i));
    end
  end
  
  if mod(count_iteration, 10) == 0
    fprintf('%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n', count_iteration, logL, abs(logL_old-logL), TOLERANCE);
  end
end

pi_i = x_rs./sum(x_rs);
pi_i = pi_i';
t = bsxfun(@rdivide, x, x_rs);

