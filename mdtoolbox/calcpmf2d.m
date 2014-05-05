function [z, xi, yi] = calcpmf2d(x, nbin)
%% calcpmf2d
% calculate 2D potential of mean force from the scattered 2D-data (using kernel density estimator)
%
%% Syntax
% [pmf, xi, yi] = calcpmf2d(data)
% [pmf, xi, yi] = calcpmf2d(data, nbin)
%
%% Description
%
%% Example
%# load p.mat;
%# nbin = 256;
%# [z, xi, yi] = calcpmf2d(p(:,[1 2]), nbin);
%# imagesc(xi,yi,z); axis xy;
%# xlabel('PC1','FontSize',45); ylabel('PC2','FontSize',45); colorbar; plot_format
%
%% See also
% calcpmf2d
%
%% References
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957. 
% 

xd = x(:,1);
yd = x(:,2);

if ~exist('nbin', 'var') || isempty(nbin)
  nbin = 256;
end

xi = linspace(min(xd),max(xd),nbin);
yi = linspace(min(yd),max(yd),nbin);
[bandwidth,z,xi,yi]=kde2d([xd yd],nbin);
xi = xi(1,:);
yi = yi(:,1)';

z(z < realmin) = NaN;
z = -log(z);
%z = -log(abs(z));
z_max = max(max(-z));
z = z + z_max;

