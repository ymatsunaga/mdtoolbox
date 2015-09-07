function [t, pi_i] = msmtransitionmatrix(c, maxiteration)
%% msmtransitionmatrix
% estimate transition probability matrix from count matrix
%
%% Syntax
%# t = msmtransitionmatrix(c);
%# t = msmtransitionmatrix(c, maxiteration);
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

c_i    = sum(c, 2);
x_i    = sum(x, 2);

if ~exist('maxiteration', 'var')
  maxiteration = 1000;
end

%% optimization by L-BFGS-B
fcn = @(x) myfunc_column(x, c, c_i, nstate);
opts.x0 = x(:);
opts.maxIts = maxiteration;
opts.maxTotalIts = 50000;
%opts.factr = 1e5;
%opts.pgtol = 1e-7;

[x, f, info] = lbfgsb(fcn, zeros(nstate*nstate, 1), Inf(nstate*nstate, 1), opts);
x = reshape(x, nstate, nstate);

x_i = sum(x, 2);
t = bsxfun(@rdivide, x, x_i);
pi_i = x_i./sum(x_i);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, g] = myfunc_column(x, c, c_i, nstate);
x = reshape(x, nstate, nstate);
[f, g] = myfunc_matrix(x, c, c_i);
g = g(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, g] = myfunc_matrix(x, c, c_i)
x_i = sum(x, 2);

% F
tmp = c .* log(bsxfun(@rdivide, x, x_i));
%index = ~(isnan(tmp));
index = (x > 0);
f = - sum(tmp(index));

% G
t = c_i./x_i;
g = (c./x) + (c'./x') - bsxfun(@plus, t, t');
index = ((x == 0) & (x' == 0));
g(index) = 0;
g = -g;

