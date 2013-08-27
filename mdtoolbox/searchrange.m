function [pair, dist] = searchrange(a, b, rcut, box)
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
%  distance 'rcut' for atoms in 'crd2'. 
%  Atom coordinates should be 3-dimensional points of vector 1 x 3natom.
%  The grid-cell algorithm by Heinz and Hünenberger (JCC, 2004)
%  is used. If the cutoff distance is large compared to the system
%  size, then the exhaustive search algorithm is used instead of
%  the grid-cell. When 'box' is given, periodic boundary condition
%  (PBC) is assumed. 
%
% * crd1 - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom' double]
% * crd2 - coordinates in 3-dimensional Cartesian space
%          [1 x 3natom'' double]
% * rcut - cutoff distance
%          [scalar double]
% * box  - box size. When this is given, PBC is assumed
%          [1 x 3 double]
% * pair - pair list. 1st column is the index of crd2 atom,
%          2nd column is the index of crd1 atom
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
%# crd = decenter(crd);
%# [pair, dist] = searchrange(crd, [0.0 0.0 0.0 3.0 4.0 5.0], 12.0);
%# 
%# % test
%# [pair2, dist2] = searchrange_exhaustive(crd, [0.0 0.0 0.0 3.0 4.0 5.0], 12.0);
%# hold off;
%# s  = sparse(pair(:, 1),  pair(:, 2),  1, 2, 164); spy(s,  'b')
%# hold on;
%# s2 = sparse(pair2(:, 1), pair2(:, 2), 1, 2, 164); spy(s2, 'r')
%# all(find(s) == find(s2))
%# all(abs(sort(dist) - sort(dist2)) < 0.00001)
%
%% See also
% searchrange_exhaustive
%
%% References
% Heinz, T.N. & Hünenberger, P.H.
% "A fast pairlist-construction algorithm for molecular simulations
% under periodic boundary conditions." 
% J Comput Chem 25, 1474–1486 (2004). 
%
%% setup variables
natom = numel(b)./3;
rcut2 = rcut^2;

%% setup cell
if nargin >= 4
  a(1:3:end) = a(1:3:end) - floor(a(1:3:end)./box(1))*box(1);
  a(2:3:end) = a(2:3:end) - floor(a(2:3:end)./box(2))*box(2);
  a(3:3:end) = a(3:3:end) - floor(a(3:3:end)./box(3))*box(3);
  b(1:3:end) = b(1:3:end) - floor(b(1:3:end)./box(1))*box(1);
  b(2:3:end) = b(2:3:end) - floor(b(2:3:end)./box(2))*box(2);
  b(3:3:end) = b(3:3:end) - floor(b(3:3:end)./box(3))*box(3);
else
  minx = min([a(1:3:end) b(1:3:end)]);
  miny = min([a(2:3:end) b(2:3:end)]);
  minz = min([a(3:3:end) b(3:3:end)]);

  a(1:3:end) = a(1:3:end) - minx;
  a(2:3:end) = a(2:3:end) - miny;
  a(3:3:end) = a(3:3:end) - minz;

  b(1:3:end) = b(1:3:end) - minx;
  b(2:3:end) = b(2:3:end) - miny;
  b(3:3:end) = b(3:3:end) - minz;

  box = zeros(1, 3);
  box(1) = max([a(1:3:end) b(1:3:end)]) + 1.0;
  box(2) = max([a(2:3:end) b(2:3:end)]) + 1.0;
  box(3) = max([a(3:3:end) b(3:3:end)]) + 1.0;
end

if any(box < (2*rcut))
  disp('calling exhaustive search...')
  [pair, dist] = rangesearch_exhaustive(crd, rcut);
  return
end

rcell = max([8.0 8.0 8.0], rcut);
ncell = floor(box./rcell);
rcell = box./ncell;

nx = floor(a(1:3:end)./rcell(1));
ny = floor(a(2:3:end)./rcell(2));
nz = floor(a(3:3:end)./rcell(3));
ma = ncell(1)*ncell(2)*nz + ncell(1)*ny + nx + 1;

nx = floor(b(1:3:end)./rcell(1));
ny = floor(b(2:3:end)./rcell(2));
nz = floor(b(3:3:end)./rcell(3));
mb = ncell(1)*ncell(2)*nz + ncell(1)*ny + nx + 1;

M = prod(ncell);

maindex  = cell(M, 1);
maindex3 = cell(M, 1);
mbindex  = cell(M, 1);
mbindex3 = cell(M, 1);
for i = 1:M
  maindex{i}  = find(ma == i);
  maindex3{i} = to3(maindex{i});
  mbindex{i}  = find(mb == i);
  mbindex3{i} = to3(mbindex{i});
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
mask_index = [-mask_index(end:-1:1) 0 mask_index];


%% calculate the distances of the atoms in masked cells
pair = zeros(natom*1000, 2);
dist = zeros(natom*1000, 1);
%pair = zeros(natom, 2);
%dist = zeros(natom, 1);
icount = 1;

for m1 = 1:M
  m1index = mbindex{m1};
  if numel(m1index) == 0
    continue
  end
  m1index3 = mbindex3{m1};

  m2 = m1 + mask_index;
  m2 = m2((1 <= m2) & (m2 <= M));
  m2index = [maindex{m2}];
  m2index3 = [maindex3{m2}];

  if nargin >= 4
    [lpair, ldist, num] = calcpair2_pbc(b(m1index3), a(m2index3), rcut, box);
  else
    [lpair, ldist, num] = calcpair2(b(m1index3), a(m2index3), rcut);
  end
  
  if num > 0
    pair(icount:(icount+num-1), :) = [m1index(lpair(:,1)')' m2index(lpair(:,2)')'];
    dist(icount:(icount+num-1)) = ldist;
    icount = icount + num;
  end
end

pair(icount:end, :) = [];
dist(icount:end) = [];

function mi_n = minimum_image(n, N)
mi_n = n - N .* sign(n) .* floor((abs(n) + 0.5*N) ./ N);


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


