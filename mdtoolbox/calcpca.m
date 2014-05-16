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
%
% * trj        - trajectory of coordinates [nstep x 3natom]
% * projection - principal components (projection of the trajectory on to principal modes) [nstep x 3natom]
% * mode       - principal modes [nstep x 3natom]
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

