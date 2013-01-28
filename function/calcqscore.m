function [qscore, dist2] = calcqscore(trj, pair, dist)
%% calcqscore
% calculate Q-score from input trajectory
%
%% Syntax
%#
%
%% Description
%
%% Example
%#
% 
%% See also
%
%% References
%

%% setup
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;
ncontact = size(pair, 1);
dist2 = zeros(nstep, ncontact);
qscore = zeros(nstep, 1);

%% calc Q-score
for istep = 1:nstep
  crd = reshape(trj(istep, :), 3, [])';
  d = arrayfun(@(x,y) sqrt(sum((crd(x, :) - crd(y, :)).^2)), pair(:,1), pair(:,2));
  qscore(istep) = sum(d < (1.15*dist)) / ncontact;
  dist2(istep, :) = d';
end

