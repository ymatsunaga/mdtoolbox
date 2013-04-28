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
% * trj            - trajectory to be clustered 
%                    [double nstep x natom3]
% * kcluster       - the number of clusters 
%                    [integer 1]
% * mass           - masses used for the calculation of rmsd 
%                    [double natom]
% * indexOfCluster - cluster index from 1 to kcluster 
%                    [integer nstep]
% * centroid       - centroids of clusters
%                    [double kcluster x natom3]
% * sumd           - sum of RMSDs from the centroid in each cluster
%                    [double nstep]
% 
%% Example
%# parm = readparm('ala.parm');
%# trj = readnetcdf('ala.nc');
%# index = find(selectid(parm.residue_id, 1:3) & ~selectname(parm.atom_name, 'H*'));
%# [indexOfCluster, centroid, sumd] = clusteringbykmeans(trj(:, to3(index)), 4, parm.mass(index));
% 
%% See also
% clusteringbykcenter
%
%% References
%
%

%% preparation
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if (nargin < 3)
  mass = [];
else
  if iscolumn(mass)
    mass = mass';
  end
end

%% remove the centers of mass of the structres
trj = decenter(trj, [], mass);

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
    rmsd = superimpose(centroid(icluster, :), trj, [], mass, [], true);
    distanceFromCentroid(:, icluster) = rmsd;
  end
  [~, indexOfCluster] = min(distanceFromCentroid, [], 2);

  % calc centroid of each cluster
  for icluster = 1:kcluster
    crd = meanstructure(trj(indexOfCluster == icluster, :), [], mass);
    centroid(icluster, :) = crd;
  end
  [~,com] = decenter(centroid, [], mass);

  % check convergence (normalized Hamming distance of indices)
  hammingDistance = double(sum(indexOfCluster ~= indexOfClusterOld))./nstep;
  fprintf(['normalized hamming distance from the previous cluster assingments: %d\n'], hammingDistance);
  if hammingDistance < 10^(-6)
    break;
  end
  
end

sumd = zeros(kcluster, 1);
for icluster = 1:kcluster
  sumd(icluster) = sum(distanceFromCentroid(indexOfCluster == icluster, icluster));
end


