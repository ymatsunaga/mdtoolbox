function [projection, mode, variances] = calcpca(trj)
%% calcpca
% perform principal component analysis of given input trajectory
%
%% Syntax
%# projection = calcpca(trj);
%# [projection, mode] = calcpca(trj);
%# [projection, mode, variances] = calcpca(trj);
%
%% Description
% This routine performs principal component analysis (PCA)
% The PCA identifies degrees of freedom which have large
% variances. 
%
% It is noted that this routine DOES NOT apply any preprocesses to
% the input trajectory, such as centering, or regularization. 
% Thus, the user may need to perform such processs before calling
% this routine. 
%
% * trj        - trajectory of coordinates [nframe x 3natom]
% * projection - principal components (projection of the trajectory on to principal modes) [nframe x 3natom]
% * mode       - principal modes [nframe x 3natom]
% * variances  - variances of principal components [3natom x 1]
%
%% Example
%# trj = readnetcdf('ak_ca.nc');
%# [~, trj] = meanstructure(trj);
%# [p, mode, variances] = calcpca(trj);
%# scatter(p(:, 1), p(:, 2), 50, 'filled');
%# xlabel('PCA 1', 'fontsize', 25);
%# ylabel('PCA 2', 'fontsize', 25);
%
%% See also
% calctica
% 

%% setup


%% covariance matrix
covar = calccovar(trj);

%% diagonalize
[eigenvector, eigenvalue] = eig(covar, 'balance');
eigenvalue = diag(eigenvalue);
[variances, index] = sort(eigenvalue, 1, 'descend');
mode = eigenvector(:, index);

%% projection
trj = bsxfun(@minus, trj, mean(trj));
projection = trj * mode;

