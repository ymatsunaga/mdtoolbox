function c = msmcountmatrix(indexOfCluster, tau, nstate)
%% msmcountmatrix
% calculate transition count matrix from a set of binned trajectory data
%
%% Syntax
%# c = msmcountmatrix(indexOfCluster);
%# c = msmcountmatrix(indexOfCluster, tau);
%# c = msmcountmatrix(indexOfCluster, tau, nstate);
%# c = msmcountmatrix(indexOfCluster, [], nstate);
%
%% Description
% calculate count matrix of transition from state i to state j during time step tau
%
%% Example
%#
% 
%% See also
%

%% setup
if ~iscell(indexOfCluster)
  indexOfCluster_noncell = indexOfCluster;
  clear indexOfCluster;
  indexOfCluster{1} = indexOfCluster_noncell;
  clear indexOfCluster_noncell;
end
ntrj = numel(indexOfCluster);

if ~exist('nstate', 'var') || isempty(nstate)
  nstate = max(cellfunc(@max, c));
  disp(sprintf('Message: nstate = %d is used.', nstate));
end

if ~exist('tau', 'var') || isempty(tau)
  tau = 1;
  disp('Message: tau = 1 is used.');
end

%% count transitions
c = sparse(nstate, nstate);

for itrj = 1:ntrj
  nframe = numel(indexOfCluster{itrj});

  index_from = 1:(nframe-tau);
  index_to   = (1+tau):nframe;
  indexOfCluster_from = indexOfCluster{itrj}(index_from);
  indexOfCluster_to   = indexOfCluster{itrj}(index_to);

  %% ignore invalid indices
  nframe = numel(indexOfCluster_from);
  s = ones(nframe, 1);

  id = (indexOfCluster_from <= 0);
  s(id) = 0;
  indexOfCluster_from(id) = 1;

  id = (indexOfCluster_to   <= 0);
  s(id) = 0;
  indexOfCluster_to(id)   = 1;

  id = isnan(indexOfCluster_from);
  s(id) = 0;
  indexOfCluster_from(id) = 1;

  id = isnan(indexOfCluster_to);
  s(id) = 0;
  indexOfCluster_to(id)   = 1;

  %% calc count matrix
  % count transitions and make count matrix C_ij by using a sparse
  % matrix
  c_itrj = sparse(indexOfCluster_from, indexOfCluster_to, s, nstate, nstate);
  c = c + c_itrj;
end

