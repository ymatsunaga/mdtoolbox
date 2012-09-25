function contactVector = calccontactvector(trj, cutoff)
%% calccontactvector
% calculate contact-map-based state-space vectors from trajectory
%
%% Syntax
%# contactVector = calccontactvector(trj);
%# contactVector = calccontactvector(trj, cutoff);
% 
%% Description
% This code creates state-space vectors from trajectory. 
% The state-space vector is caculated from the contact map 
% in which the contact of an atom pair is represented by 0 or 1. 
% This code is often used with a subsequent dimensional reduction analysis
% (e.g., principal component analysis). 
%
% * trj           - trajectory
%                   [nstep x 3natom]
% * cutoff        - threshold distance to make contact. default is 8 angstrom. 
%                   [scalar]
% * contactVector - contact map vector
%                   [nstep x natom(natom-1)/2]
%
%% Example
%# trj = readdcd('bln.dcd');
%# contactVector = calccontactvector(trj, 8.0);
%#
%# [Principal Component Analysis]
%# [pmode, p, variances, t] = princomp(contactVector);
%# scatter(p(:,1), p(:,2), 10, 'filled'); xlabel('PC1', 'FontSize', 45); ylabel('PC2', 'FontSize', 45);
%#
%# [K-means Clustering]
%# nCluster = 4;
%# [indexOfCluster, centroid] = kmeans(p(:, 1:2), nCluster, 'distance', 'sqEuclidean', 'Replicates', 5);
%# scatter(p(:, 1), p(:, 2), 5, indexOfCluster, 'filled');
%# hold on
%# scatter(centroid(:, 1), centroid(:, 2), 100, 'ko', 'filled'); xlabel('PC1 [A]', 'FontSize', 45); ylabel('PC2 [A]', 'FontSize', 45);
%# colorbar % check color-cluster correspondence
%#
%# % find the structure nearest to each centroid
%# pNearest = zeros(nCluster, 2);
%# trjNearest = zeros(nCluster, size(trj, 2));
%# for iCluster = 1:nCluster
%#   pTmp = p(indexOfCluster == iCluster, 1:2);
%#   trjTmp = trj(indexOfCluster == iCluster, :);
%#   [~, indexOfNearest] = min(sum(bsxfun(@minus, pTmp, centroid(iCluster,:)).^2, 2));
%#   pNearest(iCluster, :) = pTmp(indexOfNearest, :);
%#   trjNearest(iCluster, :) = trjTmp(indexOfNearest, :);
%# end
%# scatter(pNearest(:, 1), pNearest(:, 2), 100, 'ro', 'filled');
%# writedcd('nearest.dcd', trjNearest);
%#
%# % write the average structure of each cluster
%# trjCentroid = zeros(nCluster, size(trj, 2));
%# for i = 1:nCluster
%#   trjCentroid(i, :) = meanstructure(trj(indexOfCluster == i, :));
%# end
%# writedcd('centroid.dcd', crd);
%#
%# figure(2)
%# silhouette(contactVector, indexOfCluster, 'Hamming')
% 
%% See also
% calccontactmap, calcdistancevector
%

nstep = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

contactVector = zeros(nstep, natom*(natom-1)/2);
for istep = 1:nstep
  contactMap = calccontactmap(trj(istep, :), cutoff);
  contactVector(istep, :) = contactMap(triu(true(size(contactMap)), 1))';
end


