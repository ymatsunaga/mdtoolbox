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
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nstep' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
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
%% References
% S. Dasgupta and P. M. Long, J. Comput. Syst. Sci. 70, 555 (2005).
% J. Sun, Y. Yao, X. Huang, V. Pande, G. Carlsson, and L. J. Guibas, Learning 24, 2 (2009).
%

%% preparation
nstep = size(trj, 1);
indexOfCenter(1) = randi([1 nstep]);
indexOfCluster = ones(nstep, 1);
distPointCenter = superimpose(trj(indexOfCenter(1), :), trj, [], mass);

for i = 2:kcluster
  [~, indexOfCenter(i)] = max(distPointCenter);
  dist = superimpose(trj(indexOfCenter(i), :), trj, [], mass);
  index = dist < distPointCenter;
  if numel(index) > 0
    distPointCenter(index) = dist(index);
    indexOfCluster(index) = i;
  end
end


