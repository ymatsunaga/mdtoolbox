function [x, com] = decenter(x, index, mass)
%% decenter 
% calculate the center of mass (com) and remove the com from coordinates or velocities
%
%% Syntax
%# [trj, com] = decenter(trj);
%# [trj, com] = decenter(trj, index_atom);
%# [trj, com] = decenter(trj, index_atom, mass);
%# [trj, com] = decenter(trj, [], mass);
%
%% Description
% Calculate the center of 'mass' from given coordinates
% specified by 'index'. 
% When 'index' is omitted, all atoms are assumed. 
% When 'mass' is omitted, uniform mass weights are assumed. 
%
% * trj         - XYZ coordinates of atoms in order
%                 (x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom))
%                 [nstep x natom3 double]
% * index_atom  - index of atoms from which the center of mass are
%                 calculated [1 x n integer]
% * mass        - atom masses [1 x natom double]
% * trj(output) - XYZ coordinates of atoms where the centers of mass
%                 are removed, i.e., the com is always equal to the origin. 
%                 [nstep x natom3 double]
% * com         - centers of mass [nstep x 3]
%
%% Example
%# trj = readdcd('ak.dcd');
%
%% See also
% orient, superimpose
%

%% setup
nstep = size(x, 1);
natom3 = size(x, 2);
natom = natom3/3;
com = zeros(nstep, 3);

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

if ~exist('mass', 'var')
  mass = [];
else
  assert(isequal(natom, numel(mass)), ...
         ['sizes of coordinates and masses are not same'])
end

%% calculate the center of mass
indexx = 3.*(index-1) + 1;
indexy = 3.*(index-1) + 2;
indexz = 3.*(index-1) + 3;

if isempty(mass)

  totalMass = numel(index);
  com(:, 1) = sum(x(:, indexx), 2) ./ totalMass;
  com(:, 2) = sum(x(:, indexy), 2) ./ totalMass;
  com(:, 3) = sum(x(:, indexz), 2) ./ totalMass;

else

  if iscolumn(mass)
    mass = mass';
  end

  mass = mass(index);
  totalMass = sum(mass);
  com(:, 1) = sum(bsxfun(@times, mass, x(:, indexx)), 2) ./ totalMass;
  com(:, 2) = sum(bsxfun(@times, mass, x(:, indexy)), 2) ./ totalMass;
  com(:, 3) = sum(bsxfun(@times, mass, x(:, indexz)), 2) ./ totalMass;

end

%% subtract the center of mass
x(:, 1:3:end) = bsxfun(@minus, x(:, 1:3:end), com(:, 1));
x(:, 2:3:end) = bsxfun(@minus, x(:, 2:3:end), com(:, 2));
x(:, 3:3:end) = bsxfun(@minus, x(:, 3:3:end), com(:, 3));

