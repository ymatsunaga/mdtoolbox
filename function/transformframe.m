function [emode, variances, covar, covar_atom] = transformframe(index_fixeddomain, external_mode, covar);
%% transformframe
% transform from the Eckart frame to a non-Eckart frame
%
%% Syntax
%# covar = transformframe(index_fixeddomain, external_mode, covar)
%# [covar, covar_atom] = transformframe(index_fixeddomain, external_mode, covar);
%
%% Description
%  This code transforms the covariance matrix obtained by normal
%  mode analysis or MD simulation in the Eckart frame 
%  to the covariance matrix in a non-Eckart frame. 
%  The transformation formula in Fuchigami et al. JCP 2010 is used. 
%
% * index_fixeddomain - index of fixed domain atom.
%                       [vector]
% * external_mode     - vector of external degrees of freedom.
%                       [3natom x 6 double]
% * emode             - modes in term of new frame.
%                       [3natom x 3natom double]
% * variances         - variances in new modes in descending order.
%                       [3natom x 1 double]
% * covar             - covariance matrix of anisotropic fluctuations. 
%                       covar(i,j) = <q_i q_j>
%                       [3natom x 3natom double]
% * covar_atom        - covariance matrix of isotropic fluctuations. 
%                       covar_atom(i,j) = <x_i x_j + y_i y_j + z_i z_j>
%                       [natom x natom double]
%
%% Example
%# % frame-transformation of T4 lysozyme
%# crd = readpdb('256l_ca.pdb');
%# [emode, frequency, covar, covar_atom] = anm(crd, 8.0);
%# [emode2, variances2, covar2, covar2_atom] = transformframe([1:11 77:164], emode(:,(end-5):end), covar);
%# imagesc(covar2_atom); axis xy; colorbar; xlabel('residue','FontSize',40); ylabel('residue','FontSize',40); formatplot2
%# plot(diag(covar2_atom)); xlabel('residue','FontSize',40); ylabel('variance [a.u.]','FontSize',40); plot format
% 
%% See also
% anm, anm_sym
%
%% References
% S. Fuchigami et al., J. Chem. Phys. 132, 104109 (2010).
%

%% initial setup
external_mode = external_mode(:, (end-5):end);
natom3 = size(external_mode, 1);
natom = natom3 / 3;

P = zeros(natom3, 1);
index_fixeddomain3 = to3(index_fixeddomain);
P(index_fixeddomain3) = 1.0;
P = diag(P);

%% transformation
U = eye(natom3) ...
  - external_mode * inv(external_mode' * P * external_mode) ...
  * external_mode' * P;

covar = U * covar * U';

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

%% diagonalize new covariance matrix and get modes and variances
[u,s,v] = svd(covar);
emode = v;
variances = diag(s);

