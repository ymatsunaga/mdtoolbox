function [rdf, center] = calcrdf(index_atom1, index_atom2, trj, box, edge);
%% calcrdf
% calculate radial distribution function of two atom types specfied index_atom1 and index_atom2
%
%% Syntax
%# [rdf, center] = calcrdf(index_atom1, index_atom2, trj, [], edge);
%# [rdf, center] = calcrdf(index_atom1, index_atom2, trj, box, edge);
%
%% Description
% Calculate radial distribution function between two atom types specfied index_atom1 and index_atom2 
%
% * index_atom1 - aton index for atom type 1 [logical index or index n]
% * index_atom2 - aton index for atom type 2 [logical index or index m]
% * trj         - trajectorys [double nstep x natom3]
% * box         - box size [double nstep x 3 or 1 x 3]
% * edge        - box size [double 1 x nbin+1]
% * rdf         - radial distribution function [double 1 x nbin]
% * center      - centers of bins [double 1 x nbin]
%
%% Example
%# psf = readpsf('run.psf');
%# index1 = selectname(psf.atom_name, 'OH2');
%# index2 = selectname(psf.name, 'OH2');
%# edge = 0.0:0.1:10.0;
%# [trj, box] = readdcd('run.dcd');
%# [rdf, center] = calcrdf(index_atom1, index_atom2, trj, box, edge);
%# plot(center, rdf); formatplot
%
%% See alo
% calcangle, calcdihedral
% 

%% setup
nstep = size(trj, 1);
npair = 0;

if ~exist('box', 'var') || isempty(box)
  error('box information is necessary for radial distribution function.');
elseif (size(box, 1) == 1)
  box = repmat(box, nstep, 1);
end

if ~exist('edge', 'var')
  edge = 0:0.1:10.0;
end
rcut  = max(edge);
center = 0.5 * (edge(2:end) + edge(1:(end-1)));
nbin  = numel(center);
  
%% calculation
count  = zeros(1, nbin);
count1 = zeros(1, nbin);
if islogical(index_atom1) & islogical(index_atom2)
  npair = nnz(index_atom1)*nnz(index_atom2) - nnz(index_atom1 & index_atom2);
else
  if islogical(index_atom1)
    index_atom1 = find(index_atom1);
  end
  if islogical(index_atom2)
    index_atom2 = find(index_atom2);
  end
  npair = numel(index_atom1)*numel(index_atom2) - numel(intersect(index_atom1, index_atom2));
end
index_atom1 = to3(index_atom1);
index_atom2 = to3(index_atom2);
crd1 = zeros(1, nnz(index_atom1));
crd2 = zeros(1, nnz(index_atom2));
for istep = 1:nstep
  crd1 = trj(istep, index_atom1);
  crd2 = trj(istep, index_atom2);
  [pair, dist] = searchrange(crd1, crd2, rcut, box(istep, :));
  index_different_pair = (dist > 10.^(-6));
  count1 = histcounts(dist(index_different_pair), edge);
  count = count + count1;
end

shell_volume = (4./3) * pi * (edge(2:end).^3 - edge(1:(end-1)).^3);
s = (npair*nstep) * sum(1.0./prod(box, 2)) * shell_volume;
rdf = count./s;

