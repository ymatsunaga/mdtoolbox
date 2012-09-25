function [emode, frequency, covar, covar_atom] = anm(crd, cutoff)
%% anm
% calculate normal modes and anisotropic fluctuations by using Anisotropic (Elastic) Network Model.
%
%% Syntax
%# emode = anm(crd);
%# emode = anm(crd, cutoff);
%# [emode, frequency] = anm(crd, cutoff);
%# [emode, frequency, covar] = anm(crd, cutoff);
%# [emode, frequency, covar, covar_atom] = anm(crd, cutoff);
%
%% Description
% This code makes a hessian matrix of the Anisotropic (Elastic) 
% Network Model from the input coordinates, and diagonalize the
% matrix to obtain the normal modes (eigenvectors) and frequencies
% (sqrt(eigenvalues)).
% Also, the covariance matrix is by calculating 
% the pseudo inverse of hessian matrix. 
% kBT = 1.0 kcal/K is assumed.
%
% * crd         - coordinates of atoms. 
%                 [1 x 3natom double]
% * cutoff      - cutoff distance of the model. default is 10 angstrom. 
%                 [scalar]
% * emode       - normal modes, each column vector corresponds to a mode.
%                 1st column vector emode(:,1) is the lowest frequency mode.
%                 The external modes (translations and rotations)
%                 are emode(:,end-6:end) [3natom x 3natom double]
% * frequency   - frequencies of the normal modes, given in ascending order.
%                 frequency(end-6:end) have zero frequencies (external modes)
%                 [3natom x 1 double]
% * covar       - covariance matrix of anisotropic fluctuations. 
%                 covar(i,j) = <q_i q_j>
%                 [3natom x 3natom double]
% * covar_atom  - covariance matrix of isotropic fluctuations. 
%                 covar_atom(i,j) = <x_i x_j + y_i y_j + z_i z_j>
%                 [natom x natom double]
%
%% Example
%# pdb = readpdb('ak.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# crd = pdb.xyz(index_ca, :)';
%# crd = crd(:)';
%# [emode, frequency, covar, covar_atom] = anm(crd, 10.0);
%# plot(diag(covar_atom)); xlabel('residue','fontsize',40); ylabel('variances [a.u.]','fontsize',40); formatplot
%# imagesc(covar_atom); axis xy; xlabel('residue','fontsize',40); ylabel('residue','fontsize',40); colorbar; formatplot2
%#
%# >> crd2 = crd + emode(:, 2)'*300;
%# >> pdb.xyz = reshape(crd2, 3, [])';
%# >> writepdb('acrb2.pdb', pdb);
% 
%% See also
% anm_sym, transformframe
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

if nargin == 1
  cutoff = 10.0;
end

%% make distance matrix S
S = zeros(natom,natom);
for i = 1:natom
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);
  S(i,:) = sqrt(x_diff.^2 + y_diff.^2 + z_diff.^2);
  S(i,i) = realmax;
end

%% make kirchhoff matrix G
G = zeros(natom,natom);
for i = 1:natom
  G(i,:) = - (S(i,:) < cutoff);
  G(i,i) = - sum(G(i,:));
end

%% make Hessian matrix H from G
H = zeros(3*natom,3*natom);
for i = 1:natom
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);

  % off-diagonal block
  H(3*(i-1)+1,1:3:end) = G(i,:) .* x_diff .* x_diff ./ (S(i,:).^2);
  H(3*(i-1)+1,2:3:end) = G(i,:) .* x_diff .* y_diff ./ (S(i,:).^2);
  H(3*(i-1)+1,3:3:end) = G(i,:) .* x_diff .* z_diff ./ (S(i,:).^2);
  H(3*(i-1)+2,2:3:end) = G(i,:) .* y_diff .* y_diff ./ (S(i,:).^2);
  H(3*(i-1)+2,3:3:end) = G(i,:) .* y_diff .* z_diff ./ (S(i,:).^2);
  H(3*(i-1)+3,3:3:end) = G(i,:) .* z_diff .* z_diff ./ (S(i,:).^2);

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
[u, s, v] = svd(H);
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

