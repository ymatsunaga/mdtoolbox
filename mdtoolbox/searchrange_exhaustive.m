function [pair, dist] = searchrange_exhaustive(a, b, rcut, box)
%% searchrange
% finds all the atoms in 'crd1' that are within cutoff distance 'rcut' for atoms in 'crd2'
%
%% Syntax
%# pair = searchrange(crd1, crd2, rcut) #Non-PBC
%# [pair, dist] = calcpairlist(crd1, crd2, rcut) #Non-PBC
%# pair = calcpairlist(crd1, crd2, rcut, box) #PBC
%# [pair, dist] = calcpairlist(crd1, crd2, rcut, box) #PBC
%
%% Description
%  This code finds all the atoms in 'crd1' that are within cutoff
%  distance 'rcut' for points in 'crd2'. 
%  Input atom coordinates should be 3-dimensional points of vector 1 x 3natom.
%  When 'box' is given, periodic boundary condition (PBC) is assumed. 
%
% * crd1 - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom' double]
% * crd2 - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom'' double]
% * rcut - cutoff distance
%          [scalar double]
% * box  - box size. When this is given, PBC is assumed
%          [1 x 3 double]
% * pair - pair list. 1st column is the index of b,
%          2nd column is the index of a,
%          [n x 2 integer]
% * dist - distance of the corresponding pair
%          [n x 1 double]
%
%% Example
%# crd = readpdb('ak_ca.pdb');
%# [pair, dist] = rangesearch_exhaustive(crd, [2 3 5; 2 3 4], 8.0);
%
%% See also
% searchrange
%

%% calculate the distances of all pairs
c1 = bsxfun(@minus,b(1:3:end)',a(1:3:end));
c2 = bsxfun(@minus,b(2:3:end)',a(2:3:end));
c3 = bsxfun(@minus,b(3:3:end)',a(3:3:end));

if exist('box', 'var') && ~isempty(box)
  c1 = c1 - box(1)*round(c1./box(1));
  c2 = c2 - box(2)*round(c2./box(2));
  c3 = c3 - box(3)*round(c3./box(3));
end

dist = sqrt(c1.^2 + c2.^2 + c3.^2);

index = find(dist < rcut);

if isempty(index)
  pair = [];
  dist = [];
else
  [pair1, pair2] = ind2sub(size(dist), index);
  pair = zeros(numel(pair1), 2);
  pair(:, 1) = pair1;
  pair(:, 2) = pair2;
  dist = dist(index);
end

