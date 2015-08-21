function [indexOfCluster_out, center_out, distanceFromCenter_out] = clusterkcenters(trj, f_max, kcluster, nReplicates)
%% clusterkcenters
% K-centers clustering by using RMSD metric
%
%% Syntax
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, f_max);
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, [], kcluster);
%
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, f_max, [], nReplicates);
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, [], kcluster, nReplicates);
%
%% Description
%
% * trj             - trajectory to be clustered [nframe x natom3 double]
% * f_max           - maximum distance of samples from cluster centers
% * kcluster        - the number of clusters [scalar integer]
%
% * indexOfCluster  - cluster index from 1 to kcluster [nframe integer]
% * center          - centers of clusters [double kcluster x natom3]
% * distanceFromCenter - distance between the points and the centers of cluster [nframe double]
% 
%% Example
%# 
% 
%% See also
% clusterhybrid
%
%% References
% This function uses the method described in
% [1] S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
% [2] J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
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
    disp(sprintf('%d iteration  f_max = %f  kcluster = %d', ireplica, distanceMax, kcluster));
    if (ireplica == 1) || (distanceMax < distanceMax_out)
      distanceMax_out = distanceMax;
      indexOfCluster_out = indexOfCluster;
      center_out = trj(indexOfCenter, :);
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
    disp(sprintf('%d iteration  f_max = %f  kcluster = %d', ireplica, distanceMax, kcluster));
    if (ireplica == 1) || (kcluster < kcluster_min)
      kcluster_min = kcluster;
      indexOfCluster_out = indexOfCluster;
      center_out = trj(indexOfCenter, :);
      distanceFromCenter_out = distanceFromCenter;
    end
  end

end

%% preprocess
if exist('trj_cell', 'var')
  indexOfCluster_cell = cell(ncell, 1);
  for i = 1:ncell
    indexOfCluster_cell{i} = indexOfCluster_out((nframes(i)+1):nframes(i+1));
  end
  indexOfCluster_out = indexOfCluster_cell;

  distanceFromCenter_cell = cell(ncell, 1);
  for i = 1:ncell
    distanceFromCenter_cell{i} = distanceFromCenter_out((nframes(i)+1):nframes(i+1));
  end
  distanceFromCenter_out = distanceFromCenter_cell;
end

