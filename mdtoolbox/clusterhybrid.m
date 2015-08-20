function [indexOfCluster, indexOfCenter, distanceFromCenter] = clusterhybrid(trj, f_max)
%% clusteringbykmeans
% Hybrid clustering of K-center and K-medoid
%
%% Syntax
%# [indexOfCluster, center, sumd] = clusterhybrid(trj, kcluster, f_max)
%
%% Description
%
% * trj            - trajectory to be clustered 
%                    [double nframe x natom3]
% * kcluster       - the number of clusters 
%                    [integer scalar]
% * f_max          - maximum distance of samples from cluster centers
%                    [integer scalar]
% * indexOfCluster - cluster index from 1 to kcluster 
%                    [integer nframe]
% * center         - centers of clusters
%                    [double kcluster x natom3]
% 
%% Example
%#
% 
%% See also
% clusterkcenter
%
%% References
%
%

%% preparation
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;
mass = [];

%% remove the centers of mass of the structres
trj = decenter(trj, [], mass);

%% k-center
[indexOfCluster, indexOfCenter, distanceFromCenter] = clusterkcenters(trj, [], f_max);

%% k-medoids
kcluster = max(indexOfCluster);
center = trj(indexOfCenter, :);
f_max = max(distanceFromCenter);
f_medoids = mean(distanceFromCenter.^2);

f_cluster = zeros(kcluster, 1);
for icluster = 1:kcluster
  f_cluster(icluster) = sum(distanceFromCenter(indexOfCluster == icluster).^2);
end

count = 0;
disp(sprintf('%d iteration  f_max = %f  f_medoids = %f', count, f_max, f_medoids));

for i = 1:100;
  % propose new centers
  center_new = center;
  indexOfCenter_new = indexOfCenter;
  for icluster = 1:kcluster
    index = (indexOfCluster == icluster);
    index_find = find(index);
    id = index_find(randi(numel(index_find)));
    crd = trj(id, :);
    dist = superimpose(crd, trj(index, :), [], mass, [], true);
    if sum(dist.^2) < f_cluster(icluster)
      center_new(icluster, :) = crd;
      indexOfCenter_new(icluster) = id;
    end
  end
  
  % re-assign each data point to the closest cluster
  distanceFromCenter_new = Inf(nframe, 1);
  indexOfCluster_new = indexOfCluster;
  for icluster = 1:kcluster
    dist = superimpose(center_new(icluster, :), trj, [], mass, [], true);
    index = (dist < distanceFromCenter_new);
    if any(index)
      distanceFromCenter_new(index) = dist(index);
      indexOfCluster_new(index) = icluster;
    end
  end

  % check if new clusters yields bettter f_max
  f_max_new = max(distanceFromCenter_new)
  f_medoids_new = mean(distanceFromCenter_new.^2)
  if f_max_new < f_max
    f_max = f_max_new;
    f_medoids = mean(distanceFromCenter_new.^2);

    center = center_new;
    indexOfCenter = indexOfCenter_new;
    distanceFromCenter = distanceFromCenter_new;
    indexOfCluster_new = indexOfCluster;
  end

  count = count + 1;
  disp(sprintf('%d iteration  f_max = %f  f_medoids = %f', count, f_max, f_medoids));
end

