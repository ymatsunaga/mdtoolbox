function [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, kcluster, mass)
%% clusteringbykmeans_msd
% K-center clustering by using rmsd measure
%
%% Syntax
%# ndexOfCenter = clusteringbykcenter(trj, kcluster)
%# ndexOfCenter = clusteringbykcenter(trj, kcluster, mass)
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj, kcluster, mass)
%
%% Description
%
% * trj            - trajectory [nstep x natom3 double]
% * kcluster       - clusters [scalar integer]
% * mass           - masses used for the calculation of rmsd [natom double]
% * indexOfCluster - cluster index from 1 to kcluster [nstep integer]
% * indexOfCenter  - time index of center coordinates [kcluster integer]
% 
%% Example
%# parm = readparm('ala.parm');
%# trj = readnetcdf('ala.nc');
%# index = find(selectid(parm.residue_id, 1:3) & ~selectname(parm.atom_name, 'H*'))
%# [indexOfCluster, indexOfCenter] = clusteringbykcenter(trj(:, to3(index)), 4, parm.mass(index));
% 
%% See also
% clusteringbykmeans, clusteringbyinformation
%
%% References
% S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
% J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
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

%% clustering
% create the centroid of the 1st cluster
indexOfCenter(1) = randi([1 nstep]);
% at first, all points belong to the 1st cluster
indexOfCluster = ones(nstep, 1);
% distance between the points and the 1st centroid
distPointCenter = superimpose(trj(indexOfCenter(1), :), trj, [], mass);

for i = 2:kcluster
  [~, indexOfCenter(i)] = max(distPointCenter);
  dist = superimpose(trj(indexOfCenter(i), :), trj, [], mass);
  index = dist < distPointCenter;
  if numel(index) > 0
    % updated if the dist to a new cluster is smaller than the previous one
    distPointCenter(index) = dist(index);
    indexOfCluster(index) = i;
  end
end


