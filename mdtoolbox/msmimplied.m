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

lambda = zeros(numel(tau), 10);
for i = 1:numel(tau)
  disp(sprintf('calculating tau = %d', tau(i)));
  C = msmcountmatrix(indexOfCluster, tau(i));
  [C2, indexOfCluster2] = msmtarjan(C, indexOfCluster);
  [T, pi_i] = msmtransitionmatrix(C2, 10^(-10));
  l = eig(full(T));
  l = sort(real(l), 1, 'descend');
  lambda(i, :) = l(2:10)';
end
implied_timescale = -bsxfun(@rdivide, tau, log(lambda));

