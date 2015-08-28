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

c_i    = sum(c, 2);
x_i    = sum(x, 2);

if ~exist('TOLERANCE', 'var')
  %TOLERANCE = 10^(-14);
  TOLERANCE = 10^(-10);
end

%% optimization by self-consistent iteration
check_convergence = inf;

count_iteration = 0;
% if issparse(c)
%   index_nnz = find(c_sym);
%   [row, col] = find(c_sym);

%   while check_convergence > TOLERANCE
%     val = c_sym(index_nnz) ./ (c_i(row)./x_i(row) + c_i(col)./x_i(col));
%     x_new = sparse(row, col, val, nstate, nstate);
%     x_new_i = sum(x_new, 2);

%     count_iteration = count_iteration + 1;
%     check_convergence = norm(full(x_i./sum(x_i) - x_new_i./sum(x_new_i)));
%     x   = x_new;
%     x_i = x_new_i;
%     if mod(count_iteration, 10) == 0
%       fprintf('%d iteration  delta(pi_i) = %d  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
%     end
%   end

%   val  = x(index_nnz)./x_i(row);
%   t = sparse(row, col, val, nstate, nstate);
%   pi_i = x_i./sum(x_i);

% else

  c_ii = diag(c);
  a = bsxfun(@plus, c_i, c_i') - c_sym;
  while check_convergence > TOLERANCE
    x_old = x;

    % diagonal components
    x_ii = diag(x);
    x_ii = c_ii.*(x_i - x_ii)./(c_i - c_ii);
    x(logical(eye(size(x)))) = x_ii;
    x_i = sum(x, 2);
    
    % off-diagonal components
    b = bsxfun(@times, c_i, bsxfun(@minus, x_i', x));
    b = b + b';
    b = b - c_sym .* (bsxfun(@plus, x_i, x_i') - 2*x);
    c = - c_sym .* bsxfun(@minus, x_i, x) .* bsxfun(@minus, x_i', x);
    x = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);

    x(logical(eye(size(x)))) = x_ii;
    x = 0.5 * (x + x');
    x_i = sum(x, 2);
  
    count_iteration = count_iteration + 1;
    check_convergence = max(abs(x(:) - x_old(:)))./std(x(:));
    if mod(count_iteration, 10) == 0
      x_per_x_i = bsxfun(@rdivide, x, x_i);
      id_nnz = ~isnan(x_per_x_i(:)) & (x_per_x_i(:) > eps('single'));
      logL = sum(c(id_nnz).*log(x_per_x_i(id_nnz)));
      fprintf('%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n', count_iteration, logL, check_convergence, TOLERANCE);
    end

  end
  t = bsxfun(@rdivide, x, x_i);
  pi_i = x_i./sum(x_i);

%end
  
