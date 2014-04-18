function [emode, frequency, covar, covar_atom] = anmsparse(crd, cutoff, k)
%% anm
% calculate normal modes and anisotropic fluctuations by using Anisotropic (Elastic) Network Model.
%
%% Syntax
%# emode = anmsparse(crd);
%# emode = anmsparse(crd, cutoff);
%# [emode, frequency] = anmsparse(crd, cutoff);
%# [emode, frequency, covar] = anmsparse(crd, cutoff);
%# [emode, frequency, covar, covar_atom] = anmsparse(crd, cutoff);
%
%% Description
% This routine performs a normal mode analysis of the Anisotropic 
% (Elastic) Network Model from the input coordinates, and returns
% the normal modes (eigenvectors) and the corresponding frequencies
% (sqrt(eigenvalues)).
%
% Also, the covariance matrix of the Cartesian xyz coordinates is
% calculated from the pseudo inverse of hessian matrix. 
% In the calculation, the spring constant of 1.0 kcal/mol/A^1, 
% and kBT = 1.0 kcal/mol/K are assumed.
%
% * crd         - coordinates of atoms. 
%                 [double natom3]
% * cutoff      - cutoff distance of the model. default is 10 angstrom. 
%                 [double scalar]
% * emode       - normal modes, each column vector corresponds to a mode.
%                 1st column vector emode(:,1) is the lowest frequency mode.
%                 The external modes (translations and rotations)
%                 are emode(:,end-5:end) 
%                 [double natom3 x natom3]
% * frequency   - frequencies of the normal modes, given in ascending order.
%                 frequency(end-5:end) have zero frequencies (external modes)
%                 [double natom3]
% * covar       - covariance matrix of anisotropic fluctuations. 
%                 covar(i,j) = <q_i q_j>
%                 [double natom3 x natom3]
% * covar_atom  - covariance matrix of isotropic fluctuations. 
%                 covar_atom(i,j) = <x_i x_j + y_i y_j + z_i z_j>
%                 [double natom x natom]
%
%% Example
%# % normal mode analysis of Ca-based model of Lysozyme
%# [pdb, crd] = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb = substruct(pdb, index_ca);
%# crd = crd(to3(index_ca));
%# crd = decenter(crd);
%# [emode, frequency, covar, covar_atom] = anm(crd, 8.0);
%#
%# % rmsf and covariances
%# plot(diag(covar_atom)); xlabel('residue','fontsize',40); ylabel('variances [a.u.]','fontsize',40); formatplot
%# imagesc(covar_atom); axis xy; xlabel('residue','fontsize',40); ylabel('residue','fontsize',40); colorbar; formatplot2
%#
%# % mode structure
%# pdb.xyz = reshape(crd, 3, [])';
%# writepdb('lys_ca.pdb', pdb);
%# crd1 = crd + emode(:, 1)'*50;
%# pdb.xyz = reshape(crd1, 3, [])';
%# writepdb('lys_ca1.pdb', pdb);
%# % visualize the mode structure with PyMOL and modevectors.py
% 
%% See also
% anm_sym, transformframe, writevmdvector
%
%% References
% A.R. Atilgan, et al., Biophys. J. 80, 505-515 (2001). 
%

%% initial setup
KBT = 1.0;
x = crd(1,1:3:end);
y = crd(1,2:3:end);
z = crd(1,3:3:end);
natom = length(x);

if ~exist('cutoff', 'var')
  cutoff = 10.0;
end
cutoff2 = cutoff.^2;

if ~exist('k', 'var')
  k = 10;
end

%% make kirchhoff matrix G
G = sparse(natom,natom);
S = sparse(natom,natom);
for i = 1:natom
  for j = 1:(i-1)
    x_diff = x(i) - x(j);
    y_diff = y(i) - y(j);
    z_diff = z(i) - z(j);
    d = x_diff.^2 + y_diff.^2 + z_diff.^2;
    if d < cutoff2
      G(i,j) = - 1;
      G(j,i) = - 1;
    end
  end
end

for i = 1:natom
  G(i,i) = - sum(G(i,:));
end

%% make Hessian matrix H from G
H = sparse(3*natom,3*natom);
for i = 1:natom
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);
  d = x_diff.^2 + y_diff.^2 + z_diff.^2;
  d(i) = realmax;

  % off-diagonal block
  H(3*(i-1)+1,1:3:end) = G(i,:) .* x_diff .* x_diff ./ d;
  H(3*(i-1)+1,2:3:end) = G(i,:) .* x_diff .* y_diff ./ d;
  H(3*(i-1)+1,3:3:end) = G(i,:) .* x_diff .* z_diff ./ d;
  H(3*(i-1)+2,2:3:end) = G(i,:) .* y_diff .* y_diff ./ d;
  H(3*(i-1)+2,3:3:end) = G(i,:) .* y_diff .* z_diff ./ d;
  H(3*(i-1)+3,3:3:end) = G(i,:) .* z_diff .* z_diff ./ d;

  H(3*(i-1)+2,1:3:end) = H(3*(i-1)+1,2:3:end);
  H(3*(i-1)+3,1:3:end) = H(3*(i-1)+1,3:3:end);
  H(3*(i-1)+3,2:3:end) = H(3*(i-1)+2,3:3:end);

  % diagonal block
  H(3*(i-1)+1,3*(i-1)+1) = 0;
  H(3*(i-1)+1,3*(i-1)+2) = 0;
  H(3*(i-1)+1,3*(i-1)+3) = 0;
  H(3*(i-1)+1,3*(i-1)+1) = - sum(H(3*(i-1)+1,1:3:end));
  H(3*(i-1)+1,3*(i-1)+2) = - sum(H(3*(i-1)+1,2:3:end));
  H(3*(i-1)+1,3*(i-1)+3) = - sum(H(3*(i-1)+1,3:3:end));

  H(3*(i-1)+2,3*(i-1)+1) = 0;
  H(3*(i-1)+2,3*(i-1)+2) = 0;
  H(3*(i-1)+2,3*(i-1)+3) = 0;
  H(3*(i-1)+2,3*(i-1)+1) = - sum(H(3*(i-1)+2,1:3:end));
  H(3*(i-1)+2,3*(i-1)+2) = - sum(H(3*(i-1)+2,2:3:end));
  H(3*(i-1)+2,3*(i-1)+3) = - sum(H(3*(i-1)+2,3:3:end));

  H(3*(i-1)+3,3*(i-1)+1) = 0;
  H(3*(i-1)+3,3*(i-1)+2) = 0;
  H(3*(i-1)+3,3*(i-1)+3) = 0;
  H(3*(i-1)+3,3*(i-1)+1) = - sum(H(3*(i-1)+3,1:3:end));
  H(3*(i-1)+3,3*(i-1)+2) = - sum(H(3*(i-1)+3,2:3:end));
  H(3*(i-1)+3,3*(i-1)+3) = - sum(H(3*(i-1)+3,3:3:end));
end

%% diagonalize hessian matrix H by singular value decomposition
%[v, s, flag] = eigs(H, k+6, 0);
[~, s, v] = svds(H, k+6, 0);
emode = v(:,(end-6):-1:1);
s = diag(s);
s = sqrt(s);
frequency = s((end-6):-1:1);

%% calculate covariance (= pseudo inverse of hessian matrix)
if nargout >= 3
  covar = KBT * emode * diag(1./(frequency.^2)) * (emode');
end

if nargout >= 4
  covar_atom = zeros(natom,natom);
  for i = 1:natom
    for j = 1:natom
      covar_atom(i,j) = covar(3*(i-1)+1, 3*(j-1)+1) ...
                      + covar(3*(i-1)+2, 3*(j-1)+2) ...
                      + covar(3*(i-1)+3, 3*(j-1)+3);
    end
  end
end

%% add the 6 external degrees of freeedon in the end of output
emode = [emode v(:,end:-1:(end-5))];
frequency = [frequency; s(end:-1:(end-5))];

