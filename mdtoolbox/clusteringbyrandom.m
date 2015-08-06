function [indexOfCluster, centroid, population] = clusteringbyrandom(trj, radius)
%% clusteringbyrandom
% Lyman & Zuckerman clustering by using RMSD metric
%
%% Syntax
%# indexOfCluster = clusteringbyrandom(trj, radius)
%# [indexOfCluster, centroid] = clusteringbyrandom(trj, radius)
%# [indexOfCluster, centroid, population] = clusteringbyrandom(trj, radius)
%
%% Description
% In this clustering altogithm, cluster centers (reference
% structures) are randomly chosen from the input trajectory. 
%
% * trj             - trajectory to be clustered [nframe x natom3 double]
% * radius          - radius of clusters [scalar double]
% * indexOfCluster  - cluster index from 1 to kcluster [nframe integer]
% * centroid        - centroids of clusters [double kcluster x natom3]
% * population      - populatinos of clusters calculted from the
%                     first and second halves of trajectory 
%                     [double kcluster x 2]
% 
%% Example
%# trj = readdcd('ak_ca.dcd');
%# [indexOfCluster, centroid, population] = clusteringbyrandom(trj, 3.5);
%# semilogy(population);
%# xlabel('reference structure', 'fontsize', 25);
%# ylabel('P_i', 'fontsize', 25);
%# legend('fitst half', 'second half');
% 
%% See also
% clusteringbykcenter, clusteringbykmeans, clusteringbyinformation
%
%% References
% [1] E. Lyman and D.M. Zuckerman, Biophys. J. 91, 164 (2006).
% [2] E. Lyman and D.M. Zuckerman, J. Phys. Chem. B 11, 12876 (2007).

%% preparation
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;
indexOfCluster = zeros(nframe, 1);
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

index_firsthalf = ((1:nframe) <= nframe/2);
index_secondhalf = ~index_firsthalf;
if nargout > 2
  kcluster = max(indexOfCluster);
  population = zeros(kcluster, 2);
  for icluster = 1:kcluster
    population(icluster, 1) = sum(indexOfCluster(index_firsthalf) == icluster);
  end
  population(:, 1) = population(:, 1) ./ sum(index_firsthalf);
  for icluster = 1:kcluster
    population(icluster, 2) = sum(indexOfCluster(index_secondhalf) == icluster);
  end
  population(:, 2) = population(:, 2) ./ sum(index_secondhalf);
end

