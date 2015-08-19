function [projection, mode, lambda] = calctica(trj, lagtime)
%% calctica
% perform time-structure based Independent Component Analysis (tICA)
%
%% Syntax
%# projection = calctica(trj);
%# [projection, mode] = calctica(trj);
%# [projection, mode, lambda] = calctica(trj);
%
%% Description
% This routine performs time-structure based Independent Component
% Analysis (tICA). The tICA identifies degrees of freedom which are
% most important in the sense that their motions are very slow.
% The algorithm of the tICA is based on solving a generalized
% eigenvalue problem of time-lagges covariance matrix. The user
% needs to specify the lagtime for the calculation of the
% covariance matrix. 
%
% It is noted that this routine DOES NOT apply any preprocesses to
% the input trajectory, such as centering, or least-square fitting
% to a reference structure. The user may need to perform such
% processs before calling this routine. 
%
% * trj         - trajectory of coordinates [nframe x natom3]
% * lagtime     - lag time in the unit of frames. The default is 1.
%                 [scalar integer]
% * projection  - projections to the modes (projection of the trajectory on to tICA modes) [nframe x natom3]
% * mode        - tICA modes (normalized eigenvectors) [nframe x n]
% * lambda      - corresponding eigenvalues in descending order [n x 1]
%
%% Example
%# trj = readnetcdf('ak_ca.nc');
%# [~, trj] = meanstructure(trj);
%# [p, mode, lambda] = calctica(trj, 100);
%# scatter(p(:, 1), p(:, 2), 50, 'filled');
%# xlabel('tICA 1', 'fontsize', 25);
%# ylabel('tICA 2', 'fontsize', 25);
%
%% See also
% calcpca
% 
%% References
% [1] Y. Naritomi and S. Fuchigami, J. Chem. Phys. 134, 065101 (2011).
% [2] C. R. Schwantes and V. S. Pande, J. Chem. Theory Comput. 9, 2000 (2013).
%

%% setup
if ~exist('lagtime', 'var')
  fprintf('lagtime of 1 frame is used.\n');
  lagtime = 1;
end

%% time-lagged covariance matrix
covar0 = calccovar(trj, 0);
covar = calccovar(trj, lagtime);

%% symmetrize the time-lagged covariance matrix
covar = 0.5*(covar + covar');

%% calc pseudo-inverse of covar0
%covar0_inv = pinv(covar0);

%% remove degeneracy for solving the generalized eigenvalue problem
[pmode, variances] = eig(covar0, 'balance');
variances = diag(variances);
[variances, index] = sort(variances, 1, 'descend');
pmode = pmode(:, index);

index = variances > 10^(-6);
pmode = pmode(:, index);
covar0 = pmode'*covar0*pmode;
covar = pmode'*covar*pmode;

%% solve the generalized eigenvalue problem
[mode, lambda] = eig(covar, covar0, 'chol');
lambda = diag(lambda);
[lambda, index] = sort(lambda, 1, 'descend');
mode = mode(:, index);

%% projection
trj = bsxfun(@minus, trj, mean(trj));
mode = pmode * mode;
projection = trj * mode;

%% normalize mode vectors
fac = sqrt(sum(mode.^2));
mode = bsxfun(@rdivide, mode, fac);

