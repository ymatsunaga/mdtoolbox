function [index_A2B, index_B2A] = assigntransitionpath(data, regionA, regionB)
%% assigntransitionpath
% assign the portion of forward/backward transition paths in given 1-dimensional data
%
%% Syntax
%# [index_A2B, index_B2A] = assigntransitionpath(data, regionA, regionB)
%
%% Description
% This routine assigns the portion of forward (A->B) / backward (A<-B) 
% transition paths in given 1-dimensional data. 
%
% * data     - 1-dimensional trajectory or data set to be assigned 
%              [nstep x 1 double]
% * regionA  - upper bound for state A (reactant)
%              [scalar double]
% * regionB  - lower bound for state B (product)
%              [scalar double]
%
% * indexA2B - logical indices for forward transition path from A (reactant) to B (product)
%              if index(i) == 1, i-th step belongs to transition path
%              [nstep x 1 logical index]
% * indexB2A - logical indices for backward transition path from B (product) to A (reactant)
%              if index(i) == 1, i-th step belongs to transition path
%              [nstep x 1 logical index]
%
%% Example
%# [index_AtoB, index_BtoA] = assigntransitionpath(rmsd, 2, 8);
% 
%% See also
% assign1dbins, assign2dbins, assignvoronoi
%

%% setup
nstep = numel(data);
index_A2B = false(nstep, 1);
index_B2A = false(nstep, 1);

% index for A or B state
indexA = data < regionA;
indexB = regionB < data;
indexI = ~(indexA | indexB);

%% change point detection
% 1 for A/B -> I, -1 for I -> A/B, in 1 step evolution
index_change = indexI(2:end) - indexI(1:end-1);

% just before AB -> I
index_AB2I = find(index_change == 1);
% just before I -> AB
index_I2AB = find(index_change == -1);

%% exclude exceptional cases
% exclude the case where the initial step is already in I
if index_I2AB(1) < index_AB2I(1)
  % I2AB transition occurs before AB2I transition
  index_I2AB(1) = [];
end

% exclude the case where the final step remains in I
if index_I2AB(end) < index_AB2I(end) 
  % after AB2I transition, no I2AB transition occur
  index_AB2I(end) = [];
end

assert(numel(index_I2AB) == numel(index_AB2I), ['inward and outward crossing events do not match']);

%% check wheather I state index is transition path or not
for i = 1:numel(index_AB2I)
  istep = index_AB2I(i);
  jstep = index_I2AB(i);
  if indexA(istep) && indexB(jstep+1)
    index_A2B((istep+1):(jstep)) = true;
  elseif indexB(istep) && indexA(jstep+1)
    index_B2A((istep+1):(jstep)) = true;
  % elseif indexA(istep) && indexA(jstep+1)
  %   disp('A -> I -> A')
  % elseif indexB(istep) && indexB(jstep+1)
  %   disp('B -> I -> B')
  % else
  %   error('assignment error');
  end
end

