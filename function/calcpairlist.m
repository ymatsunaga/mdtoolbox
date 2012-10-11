function [pair, dist] = calcpairlist(crd, rcut, box)
%% calcpairlist
% make a pairlist by searching pairs within a cutoff distance
%
%% Syntax
%# pair = calcpairlist(crd, rcut) #Non-PBC
%# [pair, dist] = calcpairlist(crd, rcut) #Non-PBC
%# pair = calcpairlist(crd, rcut, box) #PBC
%# [pair, dist] = calcpairlist(crd, rcut, box) #PBC
%
%% Description
%  This code searches all the pairs within a cutoff distance.
%  Input coordinates should be 3-dimensional points of vector 1 x 3natom.
%  The grid-cell algorithm proposed by Heinz and Hünenberger (JCC, 2004)
%  is used. If the cutoff distance is large compared to the system
%  size, then the exhaustive search algorithm is used instead of
%  the grid-cell.
%
% * crd  - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom double]
% * rcut - cutoff distance
%          [scalar double]
% * box  - box size. When this is given, PBC is assumed
%          [1 x 3 double]
% * pair - pair list.
%          [npair x 2 integer]
% * dist - distance of the corresponding pair
%          [npair x 1 double]
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
%# hold off;
%# s  = sparse([pair(:,1);  pair(:,2)],  [pair(:,2);  pair(:,1)], 1, 164, 164);  spy(s,  'b')
%# hold on;
%# s2 = sparse([pair2(:,1); pair2(:,2)], [pair2(:,2); pair2(:,1)], 1, 164, 164); spy(s2, 'r')
%# all(find(s) == find(s2))
%# all(abs(sort(dist) -sort(dist2)) < 0.00001)
%
%% See also
% calcpairlist_exhaustive
%
%% References
% Heinz, T.N. & Hünenberger, P.H.
% A fast pairlist-construction algorithm for molecular simulations
% under periodic boundary conditions. 
% J Comput Chem 25, 1474–1486 (2004). 
%

%% setup variables
natom = numel(crd)./3;
rcut2 = rcut^2;

%% setup cell
if nargin >= 3
  crd(1:3:end) = crd(1:3:end) - floor(crd(1:3:end)./box(1))*box(1);
  crd(2:3:end) = crd(2:3:end) - floor(crd(2:3:end)./box(2))*box(2);
  crd(3:3:end) = crd(3:3:end) - floor(crd(3:3:end)./box(3))*box(3);
else
  minx = min(crd(1:3:end));
  miny = min(crd(2:3:end));
  minz = min(crd(3:3:end));

  crd(1:3:end) = crd(1:3:end) - minx;
  crd(2:3:end) = crd(2:3:end) - miny;
  crd(3:3:end) = crd(3:3:end) - minz;

  box = zeros(1, 3);
  box(1) = max(crd(1:3:end)) + 1.0;
  box(2) = max(crd(2:3:end)) + 1.0;
  box(3) = max(crd(3:3:end)) + 1.0;
end

if any(box < (2*rcut))
  disp('calling exhaustive search...')
  [pair, dist] = calcpairlist_exhaustive(crd, rcut);
  return
end

%rcell = [rcut rcut rcut];
rcell = max([8.0 8.0 8.0], rcut);
ncell = floor(box./rcell);
rcell = box./ncell;

nx = floor(crd(1:3:end)./rcell(1));
ny = floor(crd(2:3:end)./rcell(2));
nz = floor(crd(3:3:end)./rcell(3));

m = ncell(1)*ncell(2)*nz + ncell(1)*ny + nx + 1;
M = prod(ncell);

mindex = cell(M, 1);
mindex3 = cell(M, 1);
for i = 1:M
  mindex{i} = find(m == i);
  mindex3{i} = to3(mindex{i});
end

%% calculate cell mask
dm = 1:(M - 1);
dmx = mod(dm, ncell(1));
dmy = floor(mod(dm, ncell(1)*ncell(2)) ./ ncell(1));
dmz = floor(dm ./ (ncell(1)*ncell(2)));

dnx = abs(minimum_image(dmx, ncell(1)));

dny = zeros(size(dmy));
logical_index = (dmx == 0) | (dmy == (ncell(2)-1) & dmz == (ncell(3)-1));
dny(logical_index) = abs(minimum_image(dmy(logical_index), ncell(2)));
dny(~logical_index) = min(abs(minimum_image(dmy(~logical_index), ncell(2))), ...
                          abs(minimum_image(dmy(~logical_index) + 1, ncell(2))));

dnz = zeros(size(dmz));
logical_index = (dmz == (ncell(3)-1)) | (dmx == 0 & dmy == 0);
dnz(logical_index) = abs(minimum_image(dmz(logical_index), ncell(3)));
dnz(~logical_index) = min(abs(minimum_image(dmz(~logical_index), ncell(3))), ...
                          abs(minimum_image(dmz(~logical_index) + 1, ncell(3))));

mask = zeros(size(dm));
mask = (max(dnx, 1) - 1).^2 * rcell(1).^2 + ...
       (max(dny, 1) - 1).^2 * rcell(2).^2 + ...
       (max(dnz, 1) - 1).^2 * rcell(3).^2 <= rcut2;
%mask = true(1, M - 1); % bug check of mask
mask_index = find(mask);

%% calculate the distances of the atoms in masked cells
pair = zeros(natom*1000, 2);
dist = zeros(natom*1000, 1);
%pair = zeros(natom, 2);
%dist = zeros(natom, 1);
icount = 1;

for m1 = 1:M
  m1index = mindex{m1};
  m1index3 = mindex3{m1};
  if numel(m1index) == 0
    continue
  elseif numel(m1index) > 1
    if nargin >= 3
      [lpair, ldist, num] = calcpair_pbc(crd(m1index3), rcut, box);
    else
      [lpair, ldist, num] = calcpair(crd(m1index3), rcut);
    end
    if num > 0
      pair(icount:(icount+num-1), :) = [m1index(lpair(:,1)')' m1index(lpair(:,2)')'];
      dist(icount:(icount+num-1)) = ldist;
      icount = icount + num;
    end
  end

  m2 = m1 + mask_index;
  m2 = m2(m2 <= M);
  if isempty(m2)
    continue
  end
  m2index = [mindex{m2}];
  m2index3 = [mindex3{m2}];
  if nargin >= 3
    [lpair, ldist, num] = calcpair2_pbc(crd(m1index3), crd(m2index3), rcut, box);
  else
    [lpair, ldist, num] = calcpair2(crd(m1index3), crd(m2index3), rcut);
  end
  if num > 0
    pair(icount:(icount+num-1), :) = [m1index(lpair(:,1)')' m2index(lpair(:,2)')'];
    dist(icount:(icount+num-1)) = ldist;
    icount = icount + num;
  end
end

pair(icount:end, :) = [];
dist(icount:end) = [];

%% sort pair index
% for i = 1:size(pair, 1)
%   if pair(i, 2) > pair(i, 1)
%     tmp = pair(i, 1);
%     pair(i, 1) = pair(i, 2);
%     pair(i, 2) = tmp;
%   end
% end

function mi_n = minimum_image(n, N)
mi_n = n - N .* sign(n) .* floor((abs(n) + 0.5*N) ./ N);


function [pair, dist, num] = calcpair(crd, rcut)
c1 = bsxfun(@minus,crd(1:3:end)',crd(1:3:end));
c2 = bsxfun(@minus,crd(2:3:end)',crd(2:3:end));
c3 = bsxfun(@minus,crd(3:3:end)',crd(3:3:end));
dist = sqrt(c1.^2 + c2.^2 + c3.^2);
index = find(tril(dist < rcut, -1));
if isempty(index)
  pair = [];
  dist = [];
  num = 0;
else
  [pair1, pair2] = ind2sub(size(dist), index);
  pair = zeros(numel(pair1), 2);
  pair(:, 1) = pair1;
  pair(:, 2) = pair2;
  dist = dist(index);
  num = numel(dist);
end


function [pair, dist, num] = calcpair2(crd1, crd2, rcut)
c1 = bsxfun(@minus,crd1(1:3:end)',crd2(1:3:end));
c2 = bsxfun(@minus,crd1(2:3:end)',crd2(2:3:end));
c3 = bsxfun(@minus,crd1(3:3:end)',crd2(3:3:end));
dist = sqrt(c1.^2 + c2.^2 + c3.^2);
index = find(dist < rcut);
if isempty(index)
  pair = [];
  dist = [];
  num = 0;
else
  [pair1, pair2] = ind2sub(size(dist), index);
  pair = zeros(numel(pair1), 2);
  pair(:, 1) = pair1;
  pair(:, 2) = pair2;
  dist = dist(index);
  num = numel(dist);
end


function [pair, dist, num] = calcpair_pbc(crd, rcut, box)
c1 = bsxfun(@minus,crd(1:3:end)',crd(1:3:end));
c1 = c1 - box(1)*round(c1./box(1));
c2 = bsxfun(@minus,crd(2:3:end)',crd(2:3:end));
c2 = c2 - box(2)*round(c2./box(2));
c3 = bsxfun(@minus,crd(3:3:end)',crd(3:3:end));
c3 = c3 - box(3)*round(c3./box(3));
dist = sqrt(c1.^2 + c2.^2 + c3.^2);
index = find(tril(dist < rcut, -1));
if isempty(index)
  pair = [];
  dist = [];
  num = 0;
else
  [pair1, pair2] = ind2sub(size(dist), index);
  pair = zeros(numel(pair1), 2);
  pair(:, 1) = pair1;
  pair(:, 2) = pair2;
  dist = dist(index);
  num = numel(dist);
end


function [pair, dist, num] = calcpair2_pbc(crd1, crd2, rcut, box)
c1 = bsxfun(@minus,crd1(1:3:end)',crd2(1:3:end));
c1 = c1 - box(1)*round(c1./box(1));
c2 = bsxfun(@minus,crd1(2:3:end)',crd2(2:3:end));
c2 = c2 - box(2)*round(c2./box(2));
c3 = bsxfun(@minus,crd1(3:3:end)',crd2(3:3:end));
c3 = c3 - box(3)*round(c3./box(3));
dist = sqrt(c1.^2 + c2.^2 + c3.^2);
index = find(dist < rcut);
if isempty(index)
  pair = [];
  dist = [];
  num = 0;
else
  [pair1, pair2] = ind2sub(size(dist), index);
  pair = zeros(numel(pair1), 2);
  pair(:, 1) = pair1;
  pair(:, 2) = pair2;
  dist = dist(index);
  num = numel(dist);
end

