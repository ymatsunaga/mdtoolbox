function [pair, dist] = calcpairlist_exhaustive(crd, rcut, box)
%% calcpairlist_exhaustive
% make a pairlist by searching pairs within a cutoff distance
%
%% Syntax
%# pair = calcpairlist(crd, rcut)
%# [pair, dist] = calcpairlist(crd, rcut)
%
%% Description
%  This code searches all the pairs within a cutoff distance.
%  Since this code performs exhaustive approach (i.e., the
%  distances between all possible pairs are calculated) is used, 
%  the speed is not as fast as the grid-cell algorithgm especially
%  for big data (calcpairlist.m). 
%
% * crd  - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom double]
% * rcut - cutoff distance
%          [scalar double]
% * box  - box size. When this is given, PBC is assumed
%          [1 x 3 double]
% * pair - pair list.
%          [n x 2 integer]
% * dist - distance of the corresponding pair
%          [n x 1 double]
%
%% Example
%# % calculate Ca atom contact pairlist
%# [pdb, crd] = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# index_ca = find(index_ca);
%# crd = crd(to3(index_ca));
%# [pair, dist] = calcpairlist(crd, 8.0);
%# 
%# % test
%# [pair2, dist2] = calcpairlist_exhaustive(crd, 8.0);
%# setdiff(pair, pair2, 'rows')
%# all(abs(sort(dist) -sort(dist2)) < 0.00001)
%
%% See also
% calcpairlist
%

%% calculate the distances of all pairs
c1 = bsxfun(@minus,crd(1:3:end)',crd(1:3:end));
c2 = bsxfun(@minus,crd(2:3:end)',crd(2:3:end));
c3 = bsxfun(@minus,crd(3:3:end)',crd(3:3:end));

if nargin >= 3
  c1 = c1 - box(1)*round(c1./box(1));
  c2 = c2 - box(2)*round(c2./box(2));
  c3 = c3 - box(3)*round(c3./box(3));
end

dist = sqrt(c1.^2 + c2.^2 + c3.^2);

%index = find(tril(dist < rcut, -1));
index = find(triu(dist < rcut, 1));
[pair1, pair2] = ind2sub(size(dist), index);
pair = [pair1 pair2];
dist = dist(index);

% %% setup
% natom = numel(crd) ./ 3;
% crd = reshape(crd, 3, natom);
% rcut2 = rcut.^2;

% %% calculate the distances of all pairs
% %pair = nchoosek(natom, 2);
% pair = zeros(natom*(natom-1)/2, 2);
% [pair(:,1), pair(:,2)] = find(tril(ones(natom), -1));

% c = bsxfun(@minus, crd(:, pair(:, 1)), crd(:, pair(:, 2)));
% dist2 = sum(c.^2)';

% logical_index = (dist2 < rcut2);
% pair = pair(logical_index, :);
% dist2 = dist2(logical_index);

