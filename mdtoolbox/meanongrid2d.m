function [z, xi, yi] = meanongrid2d(x, y, data, nbin)
% average data on a regular grid
%
% [z, xi, yi] = meanongrid2d(x, y, data, nbin)
%
% nbin = 100;
% [z,xi,yi] = meanongrid2d(p(:,1), p(:,2), nbin);
% imagesc(xi,yi,z); axis xy;
% xlabel('PC1','FontSize',45); ylabel('PC2','FontSize',45); colorbar; formatplot
%

if ~exist('nbin', 'var') || isempty(nbin)
  nbin = 100;
end

%% translation of x and y to on-grid coordinates
x_min = min(x);
y_min = min(y);

x = x - x_min;
y = y - y_min;

x_max = max(x) + eps*nbin; % offset(eps*nbin) is required 
y_max = max(y) + eps*nbin; % otherwise x_index of max(x) becomes nbin + 1

x = (x ./ x_max) * nbin + 0.5;
y = (y ./ y_max) * nbin + 0.5;

x_index = round(x);
y_index = round(y);

% calculate average on grid by using sparse matrix
spmat_value = sparse(x_index, y_index, data, nbin, nbin);
spmat_count = sparse(x_index, y_index, ones(size(data)), nbin, nbin);

z = full(spmat_value./spmat_count);
z = z';
xi = ((1:nbin)./nbin) * (x_max + eps*nbin) + x_min;
yi = ((1:nbin)./nbin) * (y_max + eps*nbin) + y_min;

