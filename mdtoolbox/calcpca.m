function [p, pmode, variances] = calcpca(trj)
%% calcpca
% perform principal component analysis of given input trajectory
%
%% Syntax
%# p = calcpca(trj);
%# [p, pmode] = calcpca(trj);
%# [p, pmode, variances] = calcpca(trj);
%
%% Description
%
% * trj       - trajectory of coordinates [nstep x 3natom]
% * p         - principal components (projection of the trajectory on to principal modes) [nstep x 3natom]
% * pmode     - principal modes [nstep x 3natom]
% * variances - variances of principal components [3natom x 1]
%
%% Example
%#
%
%% See alo
%
% 

%% setup
nstep = size(trj, 1);

%% covariance matrix
trj = bsxfun(@minus, trj, mean(trj));
covmat = (trj'*trj)./(nstep-1);  % unbiased estimates of covariances
%covmat = (trj'*trj)./nstep;     % biased estimates of covariances

%% diagonalize
[eigenVector, eigenValue] = eig(covmat);
eigenValue = diag(eigenValue);
[variances, index] = sort(eigenValue, 1, 'descend');
pmode = eigenVector(:, index);

%% projection
p = trj * pmode;

