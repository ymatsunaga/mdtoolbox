function [indexOfCluster, centroid, sumd] = clusteringbykmeans_msd(trj, kcluster, mass)
%% clusteringbykmeans_msd
% K-means clustering by using rmsd measure
%
% [indexOfCluster, centroid, sumd] = clusteringbykmeans_msd(trj, kcluster, mass)
%
% input 
% trj (nstep x natom*3): trajectory
% kcluster (integer): # of clusters
% mass (1 x natom): masses used for the calculation of rmsd
%
% output
% indexOfCluster (nstep x 1): index of cluster
% centroid (kcluster x natom*3): centroids of clusters
% sumd (kcluster x 1): sum of RMSDs from the centroid in each cluster
% 
% example:
% 

%% preparation
nstep = size(trj, 1);

% create centroid by randomly drawn from input
indexOfCentroid = randperm(nstep, kcluster);
centroid = trj(indexOfCentroid, :);

%% iteration
indexOfCluster = zeros(nstep, 1);
while true

  indexOfClusterOld = indexOfCluster;
  
  % calc distance and assign cluster-index
  distanceFromCentroid = zeros(nstep, kcluster);
  for icluster = 1:kcluster
    rmsd = superimpose(centroid(icluster, :), trj);
    distanceFromCentroid(:, icluster) = rmsd;
  end
  [~, indexOfCluster] = min(distanceFromCentroid, [], 2);

  % calc centroid of each cluster
  for icluster = 1:kcluster
    crd = meanstructure(trj(indexOfCluster == icluster, :))
    centroid(icluster, :) = crd;
  end

  % check convergence (normalized Hamming distance of indices)
  hammingDistance = double(sum(indexOfCluster ~= indexOfClusterOld))./nstep
  if hammingDistance < 10^(-6)
    break;
  end
  
end

sumd = zeros(kcluster, 1);
for icluster = 1:kcluster
  sumd(icluster) = sum(distanceFromCentroid(indexOfCluster == icluster, icluster));
end


