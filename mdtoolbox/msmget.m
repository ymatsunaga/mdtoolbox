function trjOfState = msmget(state, indexOfCluster, trj)
%% msmget
% get physical quantieis from sequence of states
%
%% Syntax
%# trjOfState = msmget(state, indexOfCluster, trj)
%
%% Description
% 
%
%% Example
%#
% 
%% See also
%
%% TODO
% 
%

%% setup
nframe = numel(state);
trjOfState = zeros(nframe, 1);

if iscell(indexOfCluster)
  indexOfCluster2 = indexOfCluster;
  indexOfCluster = [];
  for i = 1:numel(indexOfCluster2)
    indexOfCluster = [indexOfCluster; indexOfCluster2{i}];
  end
  clear indexOfCluster2;
end

if iscell(trj)
  trj2 = trj;
  trj = [];
  for i = 1:numel(trj2)
    trj = [trj; trj2{i}];
  end
  clear trj2;
end

%% get physical quantieis
for iframe = 1:nframe
  index = (indexOfCluster == state(iframe));
  index = find(index);
  nz = numel(index);
  ri = randi(nz);
  trjOfState(iframe) = trj(index(ri));
end

