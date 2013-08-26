function [p, pmode, variances] = calcpca(trj)
%% calcpca
% perform principal component analysis of given input trajectory
%
%% Syntax
%# q = calcbond(trj);
%
%% Description
%
% * trj    - trajectory [nstep x 3natom]
%
%% Example
%#
%
%% See alo
% calcangle, calcdihedral
% 

%% setup
nstep = size(trj, 1);

%% covariance matrix
trj = bsxfun(@minus, trj, mean(trj));
covmat = (trj'*trj)./(nstep-1); %unbiased
%covmat = (trj'*trj)./nstep;     %biased

%% diagonalize
[eigenVector, eigenValue] = eig(covmat);
eigenValue = diag(eigenValue);
[variances, index] = sort(eigenValue, 1, 'descend');
pmode = eigenVector(:, index);

%% projection
p = trj * pmode;

