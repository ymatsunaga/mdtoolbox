function [indexOfCluster, centroid, population] = clusteringbyrandom(trj, radius)
%% clusteringbyrandom
% Lyman & Zuckerman clustering by using RMSD metric
%
%% Syntax
%# indexOfCluster = clusteringbyrandom(trj, radius)
%
%% Description
%
% * trj             - trajectory to be clustered [nstep x natom3 double]
% * radius          - radius of clusters [scalar double]
% * indexOfCluster  - cluster index from 1 to kcluster [nstep integer]
% * centroid       - centroids of clusters [double kcluster x natom3]
% 
%% Example
%#
% 
%% See also
% clusteringbykcenter, clusteringbykmeans, clusteringbyinformation
%
%% References
% E. Lyman and D.M. Zuckerman, Biophys. J. 91, 164 (2006).
% E. Lyman and D.M. Zuckerman, J. Phys. Chem. B 11, 12876 (2007).

%% preparation
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;
indexOfCluster = zeros(nstep, 1);
centroid = [];
population = [];

%% remove the centers of mass of the structres
trj = decenter(trj);

%% clustering
icount = 0;

while true
  lindex = (indexOfCluster == 0);
  index = find(lindex);
  if isempty(index)
    break;
  end
  icount = icount + 1;
  indexOfCenter(icount) = index(randi(numel(index)));
  distPointCenter = superimpose(trj(indexOfCenter(icount), :), trj(index, :), [], [], [], true);
  indexOfCluster(index(distPointCenter < radius)) = icount;

  if nargout > 1
    centroid = [centroid; trj(indexOfCenter(icount), :)];
  end
end

fprintf('%d clusters were identified.\n', icount);

if nargout > 2
  kcluster = max(indexOfCluster);
  population = zeros(kcluster, 1);
  for icluster = 1:kcluster
    population(icluster) = sum(indexOfCluster == icluster);
  end
  population = population ./ nstep;
end

