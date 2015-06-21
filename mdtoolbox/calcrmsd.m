function rmsd = calcrmsd(ref, trj, index, mass)
%% calcrmsd
% calculate root mean square deviation (RMSD) of atomic positions
%
%% Syntax
%# rmsd = calcrmsd(ref, trj);
%# rmsd = calcrmsd(ref, trj, index);
%# rmsd = calcrmsd(ref, trj, [], mass);
%# rmsd = calcrmsd(ref, trj, index, mass);
%
%% Description
% Calculate RMSD without performing any fittings
%
% * ref    - reference coordinates [1 x (natom*3) double]
% * trj    - trajectory [nstep x (natom*3) double]
% * mass   - mass [1 x natom double] or [natom x 1 double]
% * rmsd   - RMSD [nstep x 1 double]
%
%% Example
%# rmsd = calcrmsd(ref, trj);
%
%% See alo
% superimpose
% 

%% preparation
natom = numel(ref)/3;

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

%% calc rmsd
ref = ref(to3(index));
trj = trj(:, to3(index));
mass = mass(index);

[nstep, natom3] = size(trj);
natom = natom3 / 3;

x_diff = bsxfun(@minus, trj(:, 1:3:end), ref(1:3:end));
y_diff = bsxfun(@minus, trj(:, 2:3:end), ref(2:3:end));
z_diff = bsxfun(@minus, trj(:, 3:3:end), ref(3:3:end));

diff = x_diff.^2 + y_diff.^2 + z_diff.^2;
diff = bsxfun(@times, mass, diff);
diff = sum(diff, 2);
diff = diff./sum(mass);
diff = sqrt(diff);

rmsd = diff;

