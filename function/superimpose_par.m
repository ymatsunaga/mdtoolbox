function [rmsd, trj] = superimpose_par(ref, trj, index)
%% superimpose_par
% parallel implementation of 'superimpose'
%
%% See also
% superimpose
% 

%% preparation
natom3 = size(ref, 2);
natom = natom3/3;
nstep = size(trj, 1);
mass = ones(1, natom);

if nargin == 2
  index = 1:natom;
end
index3 = to3(index);

trj = decenter(trj, index, mass);
[ref, comy] = decenter(ref, index, mass);

mass = mass(index);
massxyz = repmat(mass, 3, 1);
y = reshape(ref(1, index3), 3, length(index));
rmsd = zeros(nstep, 1);

%% superimpose
trj2 = zeros(size(trj));

matlabpool
parfor istep = 1:nstep
  % calculate rotation matrix
  x = reshape(trj(istep, index3), 3, length(index));
  rmsd(istep) = 0.5 * sum(mass.*sum(x.^2 + y.^2));
  R = (massxyz.*y) * x';
  [V, D, W] = svd(R);
  D = diag(D);
  isReflection = det(V)*det(W');
  if(isReflection < 0) 
    D(3) = -D(3);
    V(:, 3) = -V(:, 3);
  end
  rmsd(istep) = rmsd(istep) - sum(D);
  U = V*W';

  % rotate molecule
  x = reshape(trj(istep, :), 3, natom);
  x = U*x;
  % translate molecule
  x(1, :) = x(1, :) + comy(1);
  x(2, :) = x(2, :) + comy(2);
  x(3, :) = x(3, :) + comy(3);
  trj2(istep, :) = reshape(x, 1, natom3);
  
end
matlabpool close

trj = trj2;
rmsd = sqrt(2.0*abs(rmsd)./sum(mass));


