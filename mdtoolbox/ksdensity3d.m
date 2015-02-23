function [f, xi, yi, zi] = ksdensity3d(data, xi, yi, zi, bandwidth, box, weight)
%% ksdensity3d
% compute 3-dimensional kernel density estimate from 3-d data
%
%% Syntax
%# f = ksdensity3d(data);
%# [f, xi, yi, zi] = ksdensity3d(data);
%# f = ksdensity3d(data, xi, yi, zi);
%# f = ksdensity3d(data, xi, yi, zi, bandwidth);
%# f = ksdensity3d(data, xi, yi, zi, bandwidth, box);
%# f = ksdensity3d(data, xi, yi, zi, bandwidth, box, weight);
%# f = ksdensity3d(data, xi, yi, zi, bandwidth, [],  weight);
%
%% Description
% Compute 3-dimensional kernel density estimate from 
% scattered 3-dimensional data.
% For the kernel function, Gaussian function is used. 
% The bandwidth can be specfied as an argument, 
% otherwise automatically determined from data and grids. 
%
% * data      - scattered 3-dimensional data [nstep x 3 double]
% * xi        - equally spaced grid in the x-axis on which
%               the density is estimated [1 x nx double]
% * yi        - equally spaced grid in the y-axis on which
%               the density is estimated [1 x ny double]
% * zi        - equally spaced grid in the z-axis on which
%               the density is estimated [1 x ny double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 3 double]
% * box       - PBC box size if data is periodic.
%               If not given or empty, data is considered to be non-periodic. 
%               [1 x 3 double]
% * weight    - weights of observables in data. 
%               By default, uniform weight is used. 
%               [nstep x 1 double]
% * f         - density estiamtes in 3-dimensional space
%               [nx x ny x nz double]
%               
%% Example
%# data = randn(1000, 3);
%# [f, xi, yi, zi] = ksdensity3d(data);
%# isosurface(xi, yi, zi, f, 0.01); grid on;
%
%% See also
% ksdensity2d
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

if ~exist('zi', 'var') || isempty(zi)
  zi = linspace(min(data(:, 3)), max(data(:, 3)), 100);
elseif iscolumn(zi)
  zi = zi';
end

nx = numel(xi);
ny = numel(yi);
nz = numel(zi);

if ~exist('bandwidth', 'var') || isempty(bandwidth)
  bandwidth = zeros(3, 1);
  
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

  sig = median(abs(zi - median(zi))) / 0.6745;
  if sig <= 0, sig = max(zi) - min(zi); end
  if sig > 0
    bandwidth(3) = sig * (1/nstep)^(1/6);
  else
    bandwidth(3) = 1;
  end

  fprintf('bandwidth in x-axis: %f\n', bandwidth(1));
  fprintf('bandwidth in y-axis: %f\n', bandwidth(2));
  fprintf('bandwidth in z-axis: %f\n', bandwidth(3));

elseif numel(bandwidth) == 1
  bandwidth = [bandwidth bandwidth bandwidth];

  fprintf('bandwidth in x-axis: %f\n', bandwidth(1));
  fprintf('bandwidth in y-axis: %f\n', bandwidth(2));
  fprintf('bandwidth in z-axis: %f\n', bandwidth(3));

end

if ~exist('box', 'var') || isempty(box)
  is_box = false;
else
  is_box = true;
end

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nstep, 1);
  weight = weight./sum(weight);
end

%% compute the kernel density estimates
if is_box
  dx = bsxfun(@minus, data(:, 1), xi);
  dx = dx - round(dx./box(1))*box(1);
  dx2 = (dx./bandwidth(1)).^2;
  dy = bsxfun(@minus, data(:, 2), yi);
  dy = dy - round(dy./box(2))*box(2);
  dy2 = (dy./bandwidth(2)).^2;
  dz = bsxfun(@minus, data(:, 3), zi);
  dz = dz - round(dz./box(3))*box(3);
  dz2 = (dz./bandwidth(3)).^2;
else
  dx2 = (bsxfun(@minus, data(:, 1), xi)./bandwidth(1)).^2;
  dy2 = (bsxfun(@minus, data(:, 2), yi)./bandwidth(2)).^2;
  dz2 = (bsxfun(@minus, data(:, 3), zi)./bandwidth(3)).^2;
end
gaussx = exp(-0.5 * dx2)./(sqrt(2*pi).*bandwidth(1));
gaussy = exp(-0.5 * dy2)./(sqrt(2*pi).*bandwidth(2));
gaussz = exp(-0.5 * dz2)./(sqrt(2*pi).*bandwidth(3));

f = zeros(nx, ny, nz);
t2 = zeros(nx, ny);
t3 = zeros(nx, ny, nz);
vec3 = zeros(1, 1, nz);
for istep = 1:nstep
  t2 = bsxfun(@times, gaussx(istep, :)', gaussy(istep, :));
  vec3(1, 1, :) = gaussz(istep, :);
  t3 = bsxfun(@times, t2, vec3);
  f = f + t3*weight(istep);
end

