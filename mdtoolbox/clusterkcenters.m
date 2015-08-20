function [indexOfCluster_out, indexOfCenter_out, distanceFromCenter_out] = clusterkcenters(trj, kcluster, f_max, nReplicates)
%% clusterkcenters
% K-center clustering by using RMSD metric
%
%% Syntax
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, kcluster);
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, [], f_max);
%
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, kcluster, [], nReplicates);
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, [], f_max, nReplicates);
%
%% Description
%
% * trj             - trajectory to be clustered [nframe x natom3 double]
% * kcluster        - the number of clusters [scalar integer]
% * f_max           - maximum distance of samples from cluster centers
% * indexOfCluster  - cluster index from 1 to kcluster [nframe integer]
% * indexOfCenter   - time index of center coordinates [kcluster integer]
% * distanceFromCenter - distance between the points and the centers of cluster [nframe double]
% 
%% Example
%# parm = readparm('ala.parm');
%# trj = readnetcdf('ala.nc');
%# index = find(selectid(parm.residue_id, 1:3) & ~selectname(parm.atom_name, 'H*'));
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj(:, to3(index)), 4, parm.mass(index));
% 
%% See also
% clusteringbykmeans, clusteringbyinformation
%
%% References
% [1] S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
% [2] J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
%

%% preparation
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;
mass = [];

if ~exist('kcluster', 'var')
  kcluster = [];
end

if ~exist('f_max', 'var')
  f_max = [];
end

if isempty(kcluster) && isempty(f_max)
  kcluster = 10;
  disp('Message: kcluster = 10 is used.');
end

if ~isempty(kcluster) && isempty(f_max)
  disp(sprintf('Message: kcluster = %d is used.', kcluster));
end

if isempty(kcluster) && ~isempty(f_max)
  disp(sprintf('Message: f_max = %f is used.', f_max));
end

if ~isempty(kcluster) && ~isempty(f_max)
  error('please do not specify kcluster and f_max together.');
end

if ~exist('nReplicates', 'var')
  nReplicates = 1;
end

%% remove the centers of mass of the structures
trj = decenter(trj, [], mass);

%% clustering
if ~isempty(kcluster)
  indexOfCenter = zeros(kcluster, 1);
  for ireplica = 1:nReplicates
    % create the center of the 1st cluster
    indexOfCenter(1) = randi([1 nframe]);
    % at first, all points belong to the 1st cluster
    indexOfCluster = ones(nframe, 1);
    % distance of the points from the 1st center
    distanceFromCenter = superimpose(trj(indexOfCenter(1), :), trj, [], mass, [], true);

    for i = 2:kcluster
      [~, indexOfCenter(i)] = max(distanceFromCenter);
      ref = trj(indexOfCenter(i), :);
      dist = superimpose(ref, trj, [], mass, [], true);
      index = (dist < distanceFromCenter);
      if any(index)
        % updated if the dist to a new cluster is smaller than the previous one
        distanceFromCenter(index) = dist(index);
        indexOfCluster(index) = i;
      end
    end

    distanceMax = max(distanceFromCenter);
    disp(sprintf('%d iteration  kcluster = %d  f_max = %f', ireplica, kcluster, distanceMax));
    if (ireplica == 1) || (distanceMax < distanceMax_out)
      distanceMax_out = distanceMax;
      indexOfCluster_out = indexOfCluster;
      indexOfCenter_out = indexOfCenter;
      distanceFromCenter_out = distanceFromCenter;
    end
  end

else
  for ireplica = 1:nReplicates
    % create the center of the 1st cluster
    indexOfCenter(1) = randi([1 nframe]);
    % at first, all points belong to the 1st cluster
    indexOfCluster = ones(nframe, 1);
    % distance of the points from the 1st center
    distanceFromCenter = superimpose(trj(indexOfCenter(1), :), trj, [], mass, [], true);

    i = 1;
    distanceMax = max(distanceFromCenter);
    while distanceMax > f_max
      i = i + 1;
      [~, indexOfCenter(i)] = max(distanceFromCenter);
      ref = trj(indexOfCenter(i), :);
      dist = superimpose(ref, trj, [], mass, [], true);
      index = dist < distanceFromCenter;
      if any(index)
        % updated if the dist to a new cluster is smaller than the previous one
        distanceFromCenter(index) = dist(index);
        indexOfCluster(index) = i;
      end
      distanceMax = max(distanceFromCenter);
    end

    kcluster = i;
    disp(sprintf('%d iteration  kcluster = %d  f_max = %f', ireplica, kcluster, distanceMax));
    if (ireplica == 1) || (kcluster < kcluster_min)
      kcluster_min = kcluster;
      indexOfCluster_out = indexOfCluster;
      indexOfCenter_out = indexOfCenter;
      distanceFromCenter_out = distanceFromCenter;
    end
  end

end

