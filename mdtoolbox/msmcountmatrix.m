function c = msmcountmatrix(indexOfCluster, nstate, tau)
%% msmcountmatrix
% calculate transition count matrix for Markov State Model (MSM)
%
%% Syntax
%# c = msmcountmatrix(indexOfCluster);
%# c = msmcountmatrix(indexOfCluster, nstate);
%# c = msmcountmatrix(indexOfCluster, nstate, tau);
%# c = msmcountmatrix(indexOfCluster, [], tau);
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
if ~exist('nstate', 'var') || isempty(nstate)
  nstate = max(indexOfCluster);
  disp(sprintf('Message: nstate = %d is used.', nstate));
end

if ~exist('tau', 'var') || isempty(tau)
  tau = 1;
  disp('Message: tau = 1 is used.');
end

nframe  = numel(indexOfCluster);

index_from = 1:(nframe-tau);
index_to   = (1+tau):nframe;
indexOfCluster_from = indexOfCluster(index_from);
indexOfCluster_to   = indexOfCluster(index_to);

%% ignore invalid indices
nframe = numel(indexOfCluster_from);
s = ones(nframe, 1);

id = (indexOfCluster_from <= 0);
s(id) = 0;
indexOfCluster_from(id) = 1;

id = (indexOfCluster_to   <= 0);
s(id) = 0;
indexOfCluster_to(id)   = 1;

id = (indexOfCluster_from == NaN);
s(id) = 0;
indexOfCluster_from(id) = 1;

id = (indexOfCluster_to   == NaN);
s(id) = 0;
indexOfCluster_to(id)   = 1;

%% calc count matrix
% count transitions and make count matrix C_ij by using a sparse matrix
c = sparse(indexOfCluster_from, indexOfCluster_to, s, nstate, nstate);

% some regularization (see the references)
%s = (c + c') > 0;
%prior = 1;
%c = c + prior*s;
% %c = 0.5*(c + c') + prior*s;

