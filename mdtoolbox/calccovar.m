function [covar, covar_atom] = calccovar(trj, lagtime)
%% calccovar
% calculate (time-lagged) covariance-matrix from input trajectory
%
%% Syntax
%# covar = calccovar(trj);
%# covar = calccovar(trj, lagtime);
%# [covar, covar_atom] = calccovar(trj, lagtime)
%
%% Description
% This routine calculates a (time-lagged) covariance matrix from
% input trajectory.
%
% * trj        - trajectory of coordinates [nstep x 3natom]
% * lagtime    - lag time in the unit of steps. The default is 0.
%                [scalar integer]
% * covar      - covariance matrix of anisotropic fluctuations. 
%                covar(i,j) = <q_i(t) q_j(t+dt)>
%                [double natom3 x natom3]
% * covar_atom - covariance matrix of atomic fluctuations. 
%                 covar_atom(i,j) = <dx_i(t) dx_j(t+dt) + dy_i(t) dy_j(t+dt) + dz_i(t) dz_j(t+dt)>
%                 [double natom x natom]
%
%% Example
%# trj = readdcd('ak.dcd');
%# [~, trj] = meanstructure(trj);
%# [~, covar_atom] = calccovar(trj);
%# imagesc(covar_atom);
%# axis xy; axis square; colorbar;
%# xlabel('residue', 'fontsize', 25);
%# ylabel('residue', 'fontsize', 25);
%# exportas('covar');
%
%% See also
% calccorr
% 

%% setup
nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if ~exist('lagtime', 'var')
  lagtime = 0;
end

%% superimpose trajectory
%fprintf('superimposing trajectory to the average structure...\n');
%[~, trj] = meanstructure(trj);

%% covariance matrix
trj = bsxfun(@minus, trj, mean(trj));
index1 = 1:(nstep-lagtime);
index2 = (1+lagtime):nstep;
covar = (trj(index1, :)'*trj(index2, :))./(nstep-lagtime-1);  % unbiased estimates of covariances
%covar = (trj'*trj)./nstep;     % biased estimates of covariances

if nargout > 1
  covar_atom = zeros(natom, natom);
  for i = 1:natom
    for j = 1:natom
      covar_atom(i,j) = covar(3*(i-1)+1, 3*(j-1)+1) ...
                      + covar(3*(i-1)+2, 3*(j-1)+2) ...
                      + covar(3*(i-1)+3, 3*(j-1)+3);
    end
  end
end

