function [indexOfCluster, center, distanceFromCenter] = clusterhybrid(trj, f_max, max_iteration)
%% clusteringbykmeans
% Hybrid clustering with K-centers and K-medoids using RMSD metric
%
%% Syntax
%# [indexOfCluster, center, sumd] = clusterhybrid(trj, f_max)
%
%% Description
%
% * trj            - trajectory to be clustered [double nframe x natom3]
% * f_max          - maximum distance of samples from cluster centers [integer scalar]
% * max_iteration  - maximum iteration of K-medoids clustering
%
% * indexOfCluster     - cluster index from 1 to kcluster [integer nframe]
% * center             - centers of clusters [double kcluster x natom3]
% * distanceFromCenter - centers of clusters [double kcluster x natom3]
% 
%% Example
%#
% 
%% See also
% clusterkcenters
%
%% References
% This function uses the method described in
% [1] K. A. Beauchamp, G. R. Bowman, T. J. Lane, L. Maibaum, I. S. Haque, and V. S. Pande, J. Chem. Theory Comput. 7, 3412 (2011).
%

%% preparation
if iscell(trj)
  trj_cell = trj;
  
  ncell = numel(trj_cell);
  natom3 = size(trj_cell{1}, 2);
  natom = natom3/3;

  nframes  = zeros(ncell+1, 1);
  for i = 1:ncell
    nframes(i+1) = nframes(i) + size(trj_cell{i}, 1);
  end

  trj = zeros(nframes(end), natom3);
  for i = 1:ncell
    trj((nframes(i)+1):nframes(i+1), :) = trj_cell{i};
  end
else
  natom3 = size(trj, 2);
  natom = natom3/3;
end

nframe = size(trj, 1);
mass = [];

if ~exist('max_iteration', 'var')
  max_iteration = Inf;
end

%% remove the centers of mass of the structres
trj = decenter(trj, [], mass);

%% k-center
disp('K-centers clustering');
[indexOfCluster, center, distanceFromCenter] = clusterkcenters(trj, f_max);

%% k-medoids
disp('K-medoids clustering');
kcluster = max(indexOfCluster);
f_max = max(distanceFromCenter);
%f_medoids = mean(distanceFromCenter.^2);
f_medoids = mean(distanceFromCenter);
f_cluster = zeros(kcluster, 1);
for icluster = 1:kcluster
  %f_cluster(icluster) = sum(distanceFromCenter(indexOfCluster == icluster).^2);
  f_cluster(icluster) = sum(distanceFromCenter(indexOfCluster == icluster));
end

count = 0;
disp(sprintf('%d iteration  f_max = %f  f_medoids = %f', count, f_max, f_medoids));

center_new = center;
indexOfCluster_new = indexOfCluster;

check_convergence = 0;
MAX_REJECTION = 10;

while (check_convergence < MAX_REJECTION) && (count < max_iteration)
  % propose new centers
  isaccepted1 = false;
  [f_max, id_max] = max(distanceFromCenter);
  %for icluster = 1:kcluster
  for icluster = randperm(kcluster, min(kcluster, 1000))
    index = (indexOfCluster == icluster);
    index_find = find(index);
    id = index_find(randi(numel(index_find)));
    crd = trj(id, :);
    dist = superimpose(crd, trj(index, :), [], mass, [], true);
    %if (sum(dist.^2) < f_cluster(icluster)) && (max(dist) < f_max)
    if (sum(dist) < f_cluster(icluster)) && (max(dist) < f_max)
      isaccepted1 = true;
      center_new(icluster, :) = crd;
    end
  end

  % re-assign each data point to the closest cluster
  if isaccepted1
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

    % check if new cluster yields bettter f_max
    f_max_new = max(distanceFromCenter_new);
    %f_medoids_new = mean(distanceFromCenter_new.^2);
    f_medoids_new = mean(distanceFromCenter_new);
    if f_max_new <= f_max
      isaccepted2 = true;

      center = center_new;
      distanceFromCenter = distanceFromCenter_new;
      indexOfCluster = indexOfCluster_new;

      f_max = f_max_new;
      f_medoids = f_medoids_new;
      f_cluster = zeros(kcluster, 1);
      for icluster = 1:kcluster
        %f_cluster(icluster) = sum(distanceFromCenter(indexOfCluster == icluster).^2);
        f_cluster(icluster) = sum(distanceFromCenter(indexOfCluster == icluster));
      end
    else
      isaccepted2 = false;
    end
  else
    isaccepted2 = false;
  end

  count = count + 1;
  if isaccepted2
    disp(sprintf('ACCEPTED at %d iteration  f_max = %f  f_medoids = %f', count, f_max, f_medoids));
    check_convergence = 0;
  else
    disp(sprintf('REJECTED at %d iteration  f_max = %f  f_medoids = %f  (proposal  f_max = %f  f_medoids = %f)', count, f_max_new, f_medoids_new));
    check_convergence = check_convergence + 1;
  end
end
disp(sprintf('Converged (proposal rejected more than %d times in succession)', MAX_REJECTION));

%% preprocess
if exist('trj_cell', 'var')
  indexOfCluster_cell = cell(ncell, 1);
  for i = 1:ncell
    indexOfCluster_cell{i} = indexOfCluster((nframes(i)+1):nframes(i+1));
  end
  indexOfCluster = indexOfCluster_cell;

  distanceFromCenter_cell = cell(ncell, 1);
  for i = 1:ncell
    distanceFromCenter_cell{i} = distanceFromCenter((nframes(i)+1):nframes(i+1));
  end
  distanceFromCenter = distanceFromCenter_cell;
end

