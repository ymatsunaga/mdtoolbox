function corr_atom = calccorr(trj)
%% calccovar
% calculate cross-correlation matrix from input trajectory
%
%% Syntax
%# covar = calccovar(trj);
%# [covar, covar_atom] = calccovar(trj);
%# [~, covar_atom] = calccovar(trj);
%
%% Description
% This routine calculates the Pearson correlation coefficients 
% between atomi and atom j, and make a cross-correlation matrix
% from input trajectory
%
% * trj       - trajectory of coordinates [nstep x 3natom]
% * corr_atom - cross correlation matrix of atomic fluctuations. 
%       covar_atom(i,j) = <dx_i dx_j + dy_i dy_j + dz_i dz_j>
%       / sqrt(<dx_i^2 + dy_i^2 + dz_i^2><dx_j^2 + dy_j^2 + dz_j^2>)
%                 [double natom x natom]
%
%% Example
%# trj = readdcd('ak.dcd');
%# [~, trj] = meanstructure(trj);
%# corr_atom = calccorr(trj);
%# imagesc(corr_atom, [-1, 1]);
%# axis xy; axis square; colorbar;
%# xlabel('residue', 'fontsize', 25);
%# ylabel('residue', 'fontsize', 25);
%# exportas('corr');
%
%% See also
% calccovar
% 

%% setup
nstep = size(trj, 1);
natom = size(trj, 2)/3;

%% calc covariance matrix
covar_atom = calccovar(trj);

%% normalize covariance matrix
corr_atom = zeros(natom, natom);
for i = 1:natom
  for j = 1:natom
    corr_atom(i,j) = covar_atom(i,j)./sqrt(covar_atom(i,i)*covar_atom(j,j));
  end
end

