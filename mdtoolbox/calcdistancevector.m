function distanceVector = calcdistancevector(trj)
%% calcdistancevector
% calculate distance-matrix-based state-space vectors from trajectory
%
%% Syntax
%# distanceVector = calcdistancevector(trj)
% 
%% Description
% This code creates state-space vectors from trajectory. 
% The state-space vector is caculated from the distance matrix. 
% This code is often used with a subsequent dimensional reduction analysis
% (e.g., principal component analysis). 
%
% * trj            - trajectory
%                    [nstep x 3natom]
% * distanceVector - distance matrix vector
%                    [nstep x natom(natom-1)/2]
%
%% Example
%# trj = readdcd('bln.dcd');
%# distanceVector = calcdistancevector(trj);
%#
%# [Principal Component Analysis]
%# [pmode, p, variances, t] = princomp(distanceVector);
%# scatter(p(:,1), p(:,2), 10, 'filled'); xlabel('PC1', 'FontSize', 45); ylabel('PC2', 'FontSize', 45);
%#
%# [K-means Clustering]
%# indexOfCluster = kmeans(distanceVector, 3, 'distance', 'sqEuclidean', 'Replicates', 5);
%# scatter(p(:, 1), p(:, 2), 10, indexOfCluster, 'filled'); xlabel('PC1', 'FontSize', 45); ylabel('PC2', 'FontSize', 45);
%# colorbar % check color-cluster correspondence
%# crd = zeros(3, size(trj, 2));
%# for i = 1:3
%#   crd(i, :) = meanstructure(trj(indexOfCluster == i, :));
%# end
%# writedcd('clusters.dcd', crd);
%# figure(2)
%# silhouette(distanceVector, indexOfCluster, 'sqEuclidean')
% 
%% See also
% calcdistancematrix, calccontactvector
%

nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

distanceVector = zeros(nstep, natom*(natom-1)/2);
for istep = 1:nstep
  distanceMatrix = calcdistancematrix(trj(istep, :));
  distanceVector(istep, :) = distanceMatrix(triu(true(size(distanceMatrix)), 1))';
end

