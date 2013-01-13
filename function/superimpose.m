function [rmsd, trj, vel, Ucell] = superimpose(ref, trj, index, mass, vel)
%% superimpose
% least-squares fitting of structures by Kabsch's method
%
%% Syntax
%# rmsd = superimpose(ref, trj);
%# rmsd = superimpose(ref, trj, index_atom);
%# rmsd = superimpose(ref, trj, index_atom, mass);
%# rmsd = superimpose(ref, trj, [], mass);
%# [rmsd, trj] = superimpose(ref, trj, index_atom, mass);
%# [rmsd, trj, vel] = superimpose(ref, trj, index_atom, mass, vel);
%# [rmsd, trj, vel] = superimpose(ref, trj, index_atom, [], vel);
%# [rmsd, trj, vel] = superimpose(ref, trj, [], [], vel);
%
%% Description
%
% * ref        - reference structure
%                [1 x natom3]
% * trj        - trajectory fitted to reference structure
%                [nstep x natom3]
% * index_atom - index of atoms to be fitted
%                [1 x n]
% * mass       - mass
%                [1 x natom]
% * vel        - velocity
%                [1 x natom3]
% * rmsd       - root mean square deviations after fitting
%                [nstep x 1]
% 
%% Example
%# trj = readnetcdf('ak.nc');
%# ref = trj(1,:);
%# [rmsd, trj] = superimpose(ref, trj);
%# plot(rmsd)
%
%% See also
% meanstructure
%
%% References
% W. Kabsch, "A solution for the best rotation to relate two sets of vectors." 
% Acta Cryst A32: 922-923 (1976)
% W. Kabsch, "A discussion of the solution for the best rotation to relate two sets of vectors." 
% Acta Cryst A34: 827-828 (1978)
% 

%% preparation
natom3 = size(ref, 2);
natom = natom3/3;
nstep = size(trj, 1);

if (nargin < 3) | (numel(index) == 0)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
  if iscolumn(index)
    index = index';
  end
end
index3 = to3(index);

if (nargin < 4) | (numel(mass) == 0)
  mass = ones(1, natom);
else
  if iscolumn(mass)
    mass = mass';
  end
end

if (nargin < 5)
  vel = [];
end

%% remove center of mass
trj = decenter(trj, index, mass);
[ref, comy] = decenter(ref, index, mass);
if numel(vel) ~= 0
  vel = decenter(vel, index, mass);
end

mass = mass(index);
massxyz = repmat(mass, 3, 1);
y = reshape(ref(1, index3), 3, length(index));
rmsd = zeros(nstep, 1);

%% superimpose
for istep = 1:nstep
  % calculate rotation matrix
  x = reshape(trj(istep, index3), 3, length(index));
  rmsd(istep) = 0.5 * sum(mass.*sum(x.^2 + y.^2));
  R = (massxyz.*y) * x';
  [V, D, W] = svd(R);
  D = diag(D);
  is_reflection = det(V)*det(W');
  if(is_reflection < 0) 
    D(3) = -D(3);
    V(:, 3) = -V(:, 3);
  end
  rmsd(istep) = rmsd(istep) - sum(D);
  U = V*W';
  if nargout >= 4
    Ucell{istep} = U;
  end

  % rotate molecule
  x = reshape(trj(istep, :), 3, natom);
  x = U*x;
  % restore the original center of mass
  x(1, :) = x(1, :) + comy(1);
  x(2, :) = x(2, :) + comy(2);
  x(3, :) = x(3, :) + comy(3);
  trj(istep, :) = reshape(x, 1, natom3);
  
  if numel(vel) ~= 0
    % rotate velocity
    v = reshape(vel(istep, :), 3, natom);
    v = U*v;
    vel(istep, :) = reshape(v, 1, natom3);
  end
end

rmsd = sqrt(2.0*abs(rmsd)./sum(mass));


