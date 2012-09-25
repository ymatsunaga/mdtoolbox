function map = calcdensity(trj, edges)
%% calcdensity
% calculate 3-D g(r) of (solvent) atoms with 3-D grids from trajectory data
%
%% Syntax
%# map = calcdensity(trj);
%# map = calcdensity(trj, edges);
%
%% Example:
% [trj, box] = readfloattrj(natom, 'md.trj', 1:natom, true);
% box = box(end, :);
% nbins = 180;
% edges{1} = linspace(0, box(1), nbins+1);
% edges{2} = linspace(0, box(2), nbins+1);
% edges{3} = linspace(0, box(3), nbins+1);
% map = density(trj(:,index3_o), edges);
% writexplormap('holo_o.xplor', map, box);
%

nstep = size(trj, 1);
natom = size(trj, 2)./3;
nbins = length(edges{1})-1;
bulk_density = (natom./(nbins.^3));

map = zeros(nbins, nbins, nbins);
crd = zeros(1, natom*3);
xyz = zeros(natom, 3);
for istep = 1:nstep
  crd = trj(istep, :);
  xyz = reshape(crd, 3, natom)';
  [map2, edges2, mid, loc] = histcn(xyz, edges{1}, edges{2}, edges{3});
  %map = map + map2;
  map = map + map2(1:nbins, 1:nbins, 1:nbins);
end
map = map./bulk_density;
map = map./nstep;

