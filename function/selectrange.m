function logical_index = selectrange(crd, index, rcut, box)
%% selectrange
% used for atom selection. Finds all the atoms within cutoff distance from indexed atoms
%
%% Syntax
%# logical_index = selectrange(crd, index, rcut); %Non-PBC
%# logical_index = selectrange(crd, index, rcut, box); %PBC
%
%% Description
% Used for atom selection. 
% Finds all the atoms within cutoff distance from indexed atoms,
% and returns logical index which has true values for the found atoms. 
% When 'box' is given, periodic boundary condition (PBC) is assumed. 
%
% * crd           - coordinates in 3-dimensional Cartesian space
%                   [1 x 3natom double]
% * index         - indices or logical indices of atoms 
%                   [1 x m char] or [1 x 3natom logical]
% * box           - box size. When this is given, PBC is assumed
%                   [1 x 3 double]
% * logical_index - logical indices which have 'trues' for the
%                   corresponding atoms
%                   [1 x natom logical]
%
%% Example
%# % read amber parm and trajectory
%# parm = readparm('ak.parm');
%# trj = readnetcdf('ak.nc');
%#
%# % choose CA atoms
%# logical_index = selectname(parm.atom_name, 'CA  ');
%# index = find(logical_index);
%#
%# % choose all the atoms within 4 Angstrom from CA atoms of residue 1-3 at step 1
%# logical_index = selectrange(trj(1, :), index(1:3), 4.0);
%# index = find(logical_index);
%

%% setup
natom3 = numel(crd);
natom = natom3./3;
index3 = to3(index);

%% searches all the atoms within cutoff distance
if nargin >= 4
  [pair, dist] = searchrange(crd, crd(index3), rcut, box);
else
  [pair, dist] = searchrange(crd, crd(index3), rcut);
end

%% get unique indices, deleting self indices
index = pair(:, 2);
index(dist < eps) = [];
index = unique(index);

logical_index = false(1, natom);
logical_index(index) = true;


