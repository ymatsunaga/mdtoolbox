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
% prior
%

%% setup
nstate = size(c, 1);

c_sym  = c + c';
x      = c_sym;

c_i    = sum(c, 2);
x_i    = sum(x, 2);

pi_i   = [];

if ~exist('TOLERANCE', 'var')
  TOLERANCE = 10^(-14);
end

%% optimization by self-consistent iteration
check_convergence = inf;

count_iteration = 0;
if issparse(c)
  index_nnz = find(c_sym);
  [row, col] = find(c_sym);

  while check_convergence > TOLERANCE
    val = c_sym(index_nnz) ./ (c_i(row)./x_i(row) + c_i(col)./x_i(col));
    x_new = sparse(row, col, val, nstate, nstate);
    x_new_i = sum(x_new, 2);
    
    count_iteration = count_iteration + 1;
    check_convergence = norm(full(x_i./sum(x_i) - x_new_i./sum(x_new_i)));
    x   = x_new;
    x_i = x_new_i;
    if mod(count_iteration, 100) == 0
      fprintf('%d iteration  delta(pi_i) = %d  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
    end
  end

  val  = x(index_nnz)./x_i(row);
  t = sparse(row, col, val, nstate, nstate);
  pi_i = x_i./sum(x_i);
  
else
  index_nnz = find(c_i > 0);

  while check_convergence > TOLERANCE
    for i = 1:numel(index_nnz)
      j = index_nnz(i);
      x1 = x_new(j, index_nnz);
      x_new(j, index_nnz) = c_sym(j, index_nnz)./(c_i(j)./x_i(j) + (c_i(index_nnz)./x_i(index_nnz))');
    end
    x_new_i = sum(x, 2);

    count_iteration = count_iteration + 1;
    check_convergence = norm(x_i./sum(x_i) - x_new_i./sum(x_new_i));
    x   = x_new;
    x_i = x_new_i;
    if mod(count_iteration, 100) == 0
      fprintf('%d iteration  delta(pi_i) = %d  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
    end
  end
  
  t = zeros(nstate, nstate);
  for i = 1:numel(index_nnz)
    j = index_nnz(i);
    t(j, index_nnz) = x(j, index_nnz)./x_i(j);
  end
  pi_i = x_i./sum(x_i);
end

