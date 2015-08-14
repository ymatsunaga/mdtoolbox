function t = msmtransitionmatrix(c)
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
x_new  = zeros(nstate, nstate);
c_i    = sum(c, 2);

%% optimization by self-consistent iteration
TOLERANCE = 10^(-14);
check_convergence = inf;

count_iteration = 0;
if issparse(c)
  index_nnz = find(c_sym);
  [row, col] = find(c_sym);

  while check_convergence > TOLERANCE
    x_i = sum(x, 2);
    val = c_sym(index_nnz) ./ (c_i(row)./x_i(row) + c_i(col)./x_i(col));
    x_new = sparse(row, col, val, nstate, nstate);

    count_iteration = count_iteration + 1;
    check_convergence = norm(x - x_new, 1);
    x = x_new;
    fprintf('%dth iteration  delta = %d  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
  end

  x_i = sum(x, 2);
  val = x(index_nnz)./x_i(row);
  t = sparse(row, col, val, nstate, nstate);
  
else
  index_nnz = find(c_i > 0);

  while check_convergence > TOLERANCE
    x_i = sum(x, 2);
    for i = 1:numel(index_nnz)
      j = index_nnz(i);
      x1 = x_new(j, index_nnz);
      x_new(j, index_nnz) = c_sym(j, index_nnz)./(c_i(j)./x_i(j) + (c_i(index_nnz)./x_i(index_nnz))');
    end

    count_iteration = count_iteration + 1;
    check_convergence = norm(x - x_new, 1);
    x = x_new;
    fprintf('%dth iteration  delta = %d  tolerance = %e\n', count_iteration, check_convergence, TOLERANCE);
  end
  
  x_i = sum(x, 2);
  t = zeros(nstate, nstate);
  for i = 1:numel(index_nnz)
    j = index_nnz(i);
    t(j, index_nnz) = x(j, index_nnz)./x_i(j);
  end
end

