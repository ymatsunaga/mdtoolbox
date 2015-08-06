function trj = orient(trj, index, mass)
%% orient 
% orient molecule using the principal axes of inertia
%
%% Syntax
% trj = orient(trj)
% trj = orient(trj, index_atom)
% trj = orient(trj, index_atom, mass)
% trj = orient(trj, [], mass)
%
%% Description
%
% * trj         - XYZ coordinates of atoms in order
%                 (x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom))
%                 [nframe x natom3 double]
% * index_atom  - index of atoms from which the center of mass are
%                 calculated [1 x n integer]
% * mass        - atom masses [1 x natom double]
% * trj(output) - XYZ coordinates of atoms in the coordinates system
%                 of the principal axes of inertia. 
%                 [nframe x natom3 double]
%
%% Example
%# [pdb, crd] = readpdb('trap.pdb');
%# crd = orient(crd);
%# pdb.xyz = reshape(crd, 3, [])';
%# writepdb('trap_orient.pdb', pdb);
%
%% See also
% decenter, superimpose
% 

%% setup
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if ~exist('index', 'var') || isempty(index)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
  if iscolumn(index)
    index = index';
  end
end

if ~exist('mass', 'var') || isempty(mass)
  mass = ones(1, natom);
else
  if iscolumn(mass)
    mass = mass';
  end
end

assert(isequal(natom, numel(mass)), ...
       'sizes of coordinates and masses are not same')

indexx = 3.*(index-1) + 1;
indexy = 3.*(index-1) + 2;
indexz = 3.*(index-1) + 3;

%% subtract by the geometric center
trj = decenter(trj, index, mass);

%% calculate the principal axis of inertia
for iframe = 1:nframe
  mass = mass(index);
  x = trj(iframe, indexx);
  y = trj(iframe, indexy);
  z = trj(iframe, indexz);
  
  I = zeros(3,3);
  
  I(1,1) = sum(mass.*(y.^2 + z.^2));
  I(2,2) = sum(mass.*(x.^2 + z.^2));
  I(3,3) = sum(mass.*(x.^2 + y.^2));
  
  I(1,2) = - sum(mass.*(x.*y));
  I(2,1) = I(1,2);
  
  I(1,3) = - sum(mass.*(x.*z));
  I(3,1) = I(1,3);
  
  I(2,3) = - sum(mass.*(y.*z));
  I(3,2) = I(2,3);
  
  [~, ~, a] = svd(I); %a is already sorted by descending order
  p_axis = a(:, end:-1:1); %z-axis has the largest inertia
  
  % check reflection
  if det(p_axis) < 0
    p_axis(:,1) = - p_axis(:,1);
  end
  
  %% project onto the principal axis of inertia
  proj = reshape(trj(iframe, :), 3, [])' * p_axis;
  trj(iframe, 1:3:end) = proj(:, 1)';
  trj(iframe, 2:3:end) = proj(:, 2)';
  trj(iframe, 3:3:end) = proj(:, 3)';
end

