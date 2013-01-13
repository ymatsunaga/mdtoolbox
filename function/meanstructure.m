function [crd, trj, vel] = meanstructure(trj, index, mass, tolerance, vel)
%% meanstructure
% calc average structure by iterative superimpose
%
%% Syntax
%# crd = meanstructure(trj);
%# crd = meanstructure(trj);
%# [crd, trj] = meanstructure(trj);
%# [crd, trj] = meanstructure(trj, index_atom);
%# [crd, trj] = meanstructure(trj, index_atom, mass);
%# [crd, trj] = meanstructure(trj, index_atom, mass, tolerance);
%# [crd, trj] = meanstructure(trj, index_atom, mass, tolerance, vel);
%# [crd, trj, vel] = meanstructure(trj, [], [], [], vel);
%
%% Description
% This routine calculates the average structure of given
% trajectory. The algorithm superimpose the trajectories to a
% plausible average structure, and update the average structrue.
% This step is repeated until some convergence in rmsd is
% achieved. 
% The total translational and rotational motions are removed in the
% output trajectory. So, this routine may be called before calling
% structure-analysis routines, such as PCA. 
%
%% Example
%# trj = readnetcdf('ak.nc');
%# [crd, trj] = meanstructure(trj);
%
%% See also
% superimpose
%

%% initialization
ref = trj(1, :);
natom3 = numel(ref);
natom = natom3/3;
rmsd = realmax;

if (nargin < 2) | (numel(index) == 0)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
end
  
if (nargin < 3) | (numel(mass) == 0)
  mass = [];
end

if (nargin < 4) | (numel(tolerance) == 0)
  tolerance = 10^(-6);
end

if nargin < 5
  vel = [];
end

%% iterative superimpose
while rmsd > tolerance
  ref_old = ref;
  [~, trj] = superimpose(ref, trj, index, mass, vel);
  ref = mean(trj);
  rmsd = superimpose(ref_old, ref, index, mass)
end

crd = ref;
[rmsd, trj, vel] = superimpose(ref, trj, index, mass, vel);


