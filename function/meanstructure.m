function [crd, trj, vel] = meanstructure(trj, index, tolerance, vel)
%% meanstructure
% calc average structure by iterative superimposing
%
%% Syntax
%# crd = meanstructure(trj);
%# crd = meanstructure(trj);
%# [crd, trj] = meanstructure(trj);
%# [crd, trj] = meanstructure(trj, index);
%# [crd, trj] = meanstructure(trj, index, tolerance);
%# [crd, trj, vel] = meanstructure(trj, index, tolerance, vel);
%# [crd, trj, vel] = meanstructure(trj, index, [], vel);
%# [crd, trj, vel] = meanstructure(trj, [], [], vel);
%
%% Description
%
%
%% Example
%
%
%% References
% W. Kabsch, "A solution for the best rotation to relate two sets of vectors." 
% Acta Cryst A32: 922-923 (1976)
% W. Kabsch, "A discussion of the solution for the best rotation to relate two sets of vectors." 
% Acta Cryst A34: 827-828 (1978)
% 

%% initialization
ref = trj(1, :);
natom = numel(ref)/3;
rmsd = realmax;

if (nargin < 2) | (numel(index) == 0)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
end
  
if (nargin < 3) | (numel(tolerance) == 0)
  tolerance = 10^(-6);
end

if nargin < 4
  vel = [];
end

%% iterative superimposing
while rmsd > tolerance
  ref_old = ref;
  [~, trj] = superimpose(ref, trj, index, [], vel);
  ref = mean(trj);
  rmsd = superimpose(ref_old, ref, index)
end

crd = ref;
[rmsd, trj, vel] = superimpose(ref, trj, index, [], vel);

