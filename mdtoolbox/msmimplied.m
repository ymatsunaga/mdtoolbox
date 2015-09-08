function [implied_timescale, lambda] = msmimplied(indexOfCluster, tau)
%% msmimplied
% calculate implied timescales
%
%% Syntax
%# implied_timescale = msmimplied(indexOfCluster, tau)
%
%% Description
% 
%
%% Example
%#
% 
%% See also
%

lambda = zeros(numel(tau), 9);
for i = 1:numel(tau)
  disp(sprintf('calculating tau = %d', tau(i)));
  C = msmcountmatrix(indexOfCluster, tau(i));
  [C2, indexOfCluster2] = msmtarjan(C, indexOfCluster);
  [T, pi_i] = msmtransitionmatrix(full(C2), 50);
  l = eig(T);
  l = sort(real(l), 1, 'descend');
  %l = eigs(T, 10, 'lr');
  lambda(i, :) = l(2:10)';
end
implied_timescale = -bsxfun(@rdivide, tau, log(lambda));

