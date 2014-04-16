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
% This routine calculates the average structure from given
% trajectory. The algorithm superimposes the trajectories to a
% plausible average structure, then updates the average structrue.
% This process is repeated until some convergence is achieved in
% the RMSD between the average structures.
% The total translational and rotational motions are removed in the
% final output trajectory which are superimposed to the converged
% average structure. So, this routine may be useful for a preprocess
% for the subsequent structure-analysis routines, such as Principal
% Component Analysis. 
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
rmsd = realmax;

if ~exist('index', 'var')
  index = [];
end
  
if ~exist('mass', 'var')
  mass = [];
end

if ~exist('tolerance', 'var') || isempty(tolerance)
  tolerance = 10^(-6);
end

if ~exist('vel', 'var')
  vel = [];
end

%% iterative superimpose
ref = decenter(ref, index, mass);
trj = decenter(trj, index, mass);
if numel(vel) ~= 0
  vel = decenter(vel, index, mass);
end

while rmsd > tolerance
  ref_old = ref;
  [~, trj, vel] = superimpose(ref, trj, index, mass, vel, true);
  ref = mean(trj);
  ref = decenter(ref, index, mass);
  rmsd = superimpose(ref_old, ref, index, mass, [], true);
  fprintf('rmsd from the previous mean structure: %f A\n', rmsd);
end

crd = ref;
[~, trj, vel] = superimpose(ref, trj, index, mass, vel);

