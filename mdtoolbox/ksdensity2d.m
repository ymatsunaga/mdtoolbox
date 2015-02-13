function [f, xi, yi] = ksdensity2d(data, xi, yi, bandwidth, box, weight)
%% ksdensity2d
% compute 2-dimensional kernel density estimate from 2-d data
%
%% Syntax
%# f = ksdensity2d(data);
%# [f, xi, yi] = ksdensity2d(data);
%# f = ksdensity2d(data, xi, yi);
%# f = ksdensity2d(data, xi, yi, bandwidth, box);
%# f = ksdensity2d(data, xi, yi, bandwidth, box, weight);
%# f = ksdensity2d(data, xi, yi, bandwidth, [],  weight);
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
% * box       - PBC box size if data is periodic.
%               If not given or empty, data is considered to be non-periodic. 
%               [1 x 2 double]
% * weight    - weights of observables in data. 
%               By default, uniform weight is used. 
%               [nstep x 1 double]
% * f         - density estiamtes in 2-dimensional space
%               [nx x ny double]
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

nx = numel(xi);
ny = numel(yi);

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
  
  fprintf('bandwidth in x-axis: %f\n', bandwidth(1));
  fprintf('bandwidth in y-axis: %f\n', bandwidth(2));

elseif numel(bandwidth) == 1
  bandwidth = [bandwidth bandwidth];

  fprintf('bandwidth in x-axis: %f\n', bandwidth(1));
  fprintf('bandwidth in y-axis: %f\n', bandwidth(2));

end

if ~exist('box', 'var') || isempty(box)
  is_box = false;
else
  is_box = true;
end

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nstep, 1);
end
weight = weight./sum(weight);

%% compute the kernel density estimates
if is_box
  dx = bsxfun(@minus, data(:, 1), xi);
  dx = dx - round(dx./box(1))*box(1);
  dx2 = (dx./bandwidth(1)).^2;
  dy = bsxfun(@minus, data(:, 2), yi);
  dy = dy - round(dy./box(2))*box(2);
  dy2 = (dy./bandwidth(2)).^2;
else
  dx2 = (bsxfun(@minus, data(:, 1), xi)./bandwidth(1)).^2;
  dy2 = (bsxfun(@minus, data(:, 2), yi)./bandwidth(2)).^2;
end
gaussx = exp(-0.5 * dx2)./(sqrt(2*pi).*bandwidth(1));
gaussy = exp(-0.5 * dy2)./(sqrt(2*pi).*bandwidth(2));

f = zeros(nx, ny);
for istep = 1:nstep
  t2 = bsxfun(@times, gaussx(istep, :)', gaussy(istep, :));
  f = f + t2*weight(istep);
end

