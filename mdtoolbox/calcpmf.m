function [z, xi] = calcpmf(x, nbin, xmin, xmax)
%% calcpmf
% calculate 1D potential of mean force from the scattered 1D-data (using kernel density estimator)
%
%% Syntax
% [pmf, xi] = calcpmf(data);
% [pmf, xi] = calcpmf(data, nbin);
% [pmf, xi] = calcpmf(data, nbin, xmin, xmax);
%
%% Description
%
%% Example
%# load p.mat
%# nbin = 256;
%# [pmf, xi] = calcpmf(x, nbin);
%# plot(xi, z,'-'); 
%# xlabel('PC1','FontSize',25); ylabel('\Delta {\itPMF}({\itk_{B}T})','FontSize',25); 
%# formatplot;
%
%% See also
% calcpmf2d
%
%% References
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957. 
% 

%% setup
%x = x(:, 1);

if ~exist('nbin', 'var') || isempty(nbin)
  nbin = 256;
end

if ~exist('xmin', 'var') || isempty(xmin)
  xmin = min(x);
end

if ~exist('xmax', 'var') || isempty(xmin)
  xmax = max(x);
end

%% calculation
%z = ksdensity(x, xi);
[bandwidth, z, xi] = kde(x, nbin);
z(z < realmin) = NaN;
z = -log(z);
z_min = min(z(:));
z = z - z_min;

