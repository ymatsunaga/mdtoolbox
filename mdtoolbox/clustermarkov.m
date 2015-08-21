function [indexOfCluster_out] = clustermarkov(t, indexOfState, r)
%% clusteringbymarkob
% Markov Cluster algorithm (MCL)
%
%% Syntax
%# [indexOfCluster] = clusteringbymarkov(t, indexOfState)
%
%% Description
% 
%% Example
% 
%% See also
% clusteringbykmeans, clusteringbyinformation
%
%% References
%

%% preparation
nstate = size(t, 1);
%nframe = size(indexOfState, 1);

if ~exist('r', 'var')
  r = 2;
end

%% clustering
t_old = t;
while true
  % expansion
  t = t^2;
  % inflation
  t = t.^r;
  t = bsxfun(@rdivide, t, sum(t));
  %t = bsxfun(@rdivide, t, sum(t, 2));
  % check convergence
  if max(max(abs(t_old - t))) < 10^(-10)
    break;
  else
    t_old = t;
  end
end

%% identify clusters
state2cluster = zeros(1, nstate);
icluster = 0;
while sum(state2cluster == 0) > 0
  icluster = icluster + 1;
  icolumn = min(find(state2cluster == 0));
  state2cluster(icolumn) = icluster;
  for i = find(state2cluster == 0)
    if max(abs(t(:, i) - t(:, icolumn))) < 10^(-10)
      state2cluster(i) = icluster;
    end
  end
end
fprintf('%d micro-states were clusterized to %d states\n', nstate, icluster);

%% cluster the states
for i = 1:icluster
  state = find(state2cluster == i);
  for j = state
    indexOfState(indexOfState == j) = i;
  end
end

indexOfCluster_out = indexOfState;

