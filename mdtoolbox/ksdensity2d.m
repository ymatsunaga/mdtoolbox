function [f, xi, yi] = ksdensity2d(data, xi, yi, bandwidth, weight)
%% ksdensity2d
% compute 2-dimensional kernel density estimate from 2-d data
%
%% Syntax
%# f = ksdensity2d(data);
%# [f, xi, yi] = ksdensity2d(data);
%# f = ksdensity2d(data, xi, yi);
%# f = ksdensity2d(data, xi, yi, bandwidth);
%# f = ksdensity2d(data, xi, yi, bandwidth, weight);
%
%% Description
% Compute 2-dimensional kernel density estimate from 
% scattered 2-dimensional data.
% For the kernel function, Gaussian function is used. 
% The bandwidth can be specfied as an argument, 
% otherwise automatically determined from data and grids. 
%
% * data      - scattered 2-dimensional data [nstep x 2 double]
% * xi        - equally spaced grid in the x-axis on which
%               the density is estimated [1 x nx double]
% * yi        - equally spaced grid in the y-axis on which
%               the density is estimated [1 x ny double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 2 double]
% * weight    - weights of observables in data. 
%               By default, uniform weight is used. 
%               [nstep x 1 double]
%               
%% Example
%# data = randn(1000, 2);
%# [f, xi, yi] = ksdensity2d(data);
%# surf(xi, yi, f);
%
%% See also
% ksdensity3d
%

%% setup
nstep = size(data, 1);

if ~exist('xi', 'var') || isempty(xi)
  xi = linspace(min(data(:, 1)), max(data(:, 1)), 100);
elseif iscolumn(xi)
  xi = xi';
end

if ~exist('yi', 'var') || isempty(yi)
  yi = linspace(min(data(:, 2)), max(data(:, 2)), 100);
elseif iscolumn(yi)
  yi = yi';
end

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nstep, 1);
end
weight = weight./sum(weight);

if ~exist('bandwidth', 'var') || isempty(bandwidth)
  bandwidth = zeros(2, 1);
  
  sig = median(abs(xi - median(xi)))/0.6745;
  if sig <= 0, sig = max(xi) - min(xi); end
  if sig > 0
    bandwidth(1) = sig * (1/nstep)^(1/6);
  else
    bandwidth(1) = 1;
  end

  sig = median(abs(yi - median(yi))) / 0.6745;
  if sig <= 0, sig = max(yi) - min(yi); end
  if sig > 0
    bandwidth(2) = sig * (1/nstep)^(1/6);
  else
    bandwidth(2) = 1;
  end
end

fprintf('bandwidth in x-axis: %f\n', bandwidth(1));
fprintf('bandwidth in y-axis: %f\n', bandwidth(2));

nx = numel(xi);
ny = numel(yi);

%% compute the kernel density estimates
dx2 = (bsxfun(@minus, data(:, 1), xi)./bandwidth(1)).^2;
dy2 = (bsxfun(@minus, data(:, 2), yi)./bandwidth(2)).^2;
gaussx = exp(-0.5 * dx2)./(sqrt(2*pi).*bandwidth(1));
gaussy = exp(-0.5 * dy2)./(sqrt(2*pi).*bandwidth(2));

f = zeros(nx, ny);
for istep = 1:nstep
  t2 = bsxfun(@times, gaussx(istep, :)', gaussy(istep, :));
  f = f + t2*weight(istep);
end

