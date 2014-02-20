function [covar, covar_atom] = calccovar(trj)
%% calccovar
% calculate covariance-matrix from trajectory
%
%% Syntax
%# covar = calccovar(trj);
%
%% Description
%
% * trj        - trajectory of coordinates [nstep x 3natom]
% * covar      - covariance matrix of anisotropic fluctuations. 
%                 covar(i,j) = <q_i q_j>
%                 [double natom3 x natom3]
% * covar_atom - covariance matrix of isotropic fluctuations. 
%                 covar_atom(i,j) = <x_i x_j + y_i y_j + z_i z_j>
%                 [double natom x natom]
%
%% Example
%#
%
%% See alo
% calcpca
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

if nargout >= 2
  covar_atom = zeros(natom,natom);
  for i = 1:natom
    for j = 1:natom
      covar_atom(i,j) = covar(3*(i-1)+1, 3*(j-1)+1) ...
                      + covar(3*(i-1)+2, 3*(j-1)+2) ...
                      + covar(3*(i-1)+3, 3*(j-1)+3);
    end
  end
end

