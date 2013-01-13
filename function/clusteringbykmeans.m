function [indexOfCluster, centroid, sumd] = clusteringbykmeans(trj, kcluster, mass)
%% clusteringbykmeans
% K-means clustering by using rmsd measure
%
%% Syntax
%# indexOfCluster = clusteringbykmeans(trj, kcluster)
%# indexOfCluster = clusteringbykmeans(trj, kcluster, mass)
%# [indexOfCluster, centroid] = clusteringbykmeans(trj, kcluster, mass)
%# [indexOfCluster, centroid, sumd] = clusteringbykmeans(trj, kcluster, mass)
%
%% Description
%
% * trj            - trajectory [nstep x natom3 double]
% * kcluster       - clusters [scalar integer]
% * mass           - masses used for the calculation of rmsd [natom double]
% * indexOfCluster - cluster index from 1 to kcluster [nstep integer]
% * centroid       - centroids of clusters [[kcluster x natom3 double]
% * sumd           - sum of RMSDs from the centroid in each cluster
% 
%% Example
%# parm = readparm('ala.parm');
%# trj = readnetcdf('ala.nc');
%# index = find(selectid(parm.residue_id, 1:3) & ~selectname(parm.atom_name, 'H*'))
%# [indexOfCluster, centroid, sumd] = clusteringbykmeans(trj(:, to3(index)), 4, parm.mass(index));
% 
%% References
% 

%% preparation
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if (nargin < 3) | (numel(mass) == 0)
  mass = ones(1, natom);
else
  if iscolumn(mass)
    mass = mass';
  end
end

%% clustering by iteration
% create centroid by randomly drawn from input
indexOfCentroid = randperm(nstep, kcluster);
centroid = trj(indexOfCentroid, :);

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
    crd = meanstructure(trj(indexOfCluster == icluster, :));
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


