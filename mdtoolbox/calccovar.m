function covar_atom = calccovar(trj)
%% calccovar
% calculate covariance-matrix from input trajectory
%
%% Syntax
%# covar_atom = calccovar(trj);
%
%% Description
% This routine calculates a covariance matrix from input
% trajectory. Prior to the calculation of the covariance,
% superimposing is automatically applied to eliminate external
% translations and rotations. In the current implementation,
% mass-weight calculation is not supported. Just assuming uniform
% mass for all atoms.
%
% * trj        - trajectory of coordinates [nstep x 3natom]
% * covar_atom - covariance matrix of atomic fluctuations. 
%                 covar_atom(i,j) = <dx_i dx_j + dy_i dy_j + dz_i dz_j>
%                 [double natom x natom]
%
%% Example
%# trj = readdcd('ak.dcd');
%# covar_atom = calccovar(trj);
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
natom = size(trj, 2)/3;

%% superimpose trajectory
fprintf('superimposing trajectory to the average structure...\n');
[~, trj] = meanstructure(trj);

%% covariance matrix
trj = bsxfun(@minus, trj, mean(trj));
covar = (trj'*trj)./(nstep-1);  % unbiased estimates of covariances
%covar = (trj'*trj)./nstep;     % biased estimates of covariances

covar_atom = zeros(natom, natom);
for i = 1:natom
  for j = 1:natom
    covar_atom(i,j) = covar(3*(i-1)+1, 3*(j-1)+1) ...
                    + covar(3*(i-1)+2, 3*(j-1)+2) ...
                    + covar(3*(i-1)+3, 3*(j-1)+3);
  end
end

