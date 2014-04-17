function ene = calcstrainenergy(ref, trj, rcut)
%% calcstrainenergy
% calculate strain energy
%
%% Syntax
%# ene = calcstrainenergy(ref, trj)
% 
%% Description
% calculate strain energy
%
% * ref   - refeference structure for elastic network model
%           [1 x 3natom]
% * trj   - elastic energies are calculated for trj
%           [nstep x 3natom]
% * ene   - elastic energies
%           [nstep x natom]
%
%% Example
%# 
% 
%% See also
% calcdistancematrix, calcdistancevector
%

%% initialization
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if ~exist('rcut', 'var') || isempty(rcut)
  rcut = 8.0;
end

%% define contacts
[pair, d0] = calcpairlist(ref, rcut);

%% calculate restrain energies
ene = zeros(nstep, natom);
d = calcbond(trj, pair);
d = bsxfun(@minus, d, d0');
d = d.^2;

for iatom = 1:natom
  index = (pair(:, 1) == iatom) | (pair(:, 2) == iatom);
  ene(:, iatom) = sum(d(:, index), 2);
end

