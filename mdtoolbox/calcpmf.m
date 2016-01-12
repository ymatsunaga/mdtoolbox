function [pmf, xi] = calcpmf(data, xi, bandwidth, box, weight)
%% calcpmf
% calculate 1D potential of mean force from the 1D-data by using kernel density estimator
%
%% Syntax
% [pmf, xi] = calcpmf(data)
% pmf = calcpmf(data, xi)
% pmf = calcpmf(data, xi, bandwidth, box, weight)
%
%% Description
%
% * data      - one-dimensional data
%               [nframe x 1 double array]
% * xi        - equally spaced grid
%               [n x 1 or 1 x n double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 1 double]
% * box       - PBC box size if data is periodic.
%               If not given or empty, data is considered to be non-periodic. 
%               [1 x 1 double]
% * weight    - weights of observables in data. 
%               By default, uniform weight is used. 
%               [nframe x 1 double]
% * pmf       - potential of mean force
%               [m x n double array]
%
%% Example
%# xi = -180:2:180; % grids
%# pmf = calcpmf(phi, xi, 3.0, 360);
%# s   = getconstants(); % get Boltzmann constant in kcal/mol/K
%# T   = 300.0;          % set temperature
%# pmf = s.KB*T*pmf;     % convert unit from KBT to kcal/mol
%# plot(xi, pmf); xlable('Phi [degrees]'); ylable('PMF [kcal/mol]');
%
%% See also
% calcpmf2d
%

%% setup
if ~exist('xi', 'var')
  xi = [];
end

if ~exist('bandwidth', 'var') || isempty(bandwidth)
  bandwidth = [];
end

if ~exist('box', 'var')
  box = [];
end

if ~exist('weight', 'var')
  weight = [];
elseif isrow(weight)
  weight = weight';
end

%% calculation
pdf = ksdensity1d(data, xi, bandwidth, box, weight);
pdf(pdf < realmin('single')) = NaN;
pmf = -log(pdf);
pmf = pmf - min(pmf(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, xi] = ksdensity1d(data, xi, bandwidth, box, weight)
% subfunction for one-dimensional kernel density estimator
%

%% setup
nframe = size(data, 1);

if ~exist('xi', 'var') || isempty(xi)
  xi = linspace(min(data(:)), max(data(:)), 100);
elseif iscolumn(xi)
  xi = xi';
end

nx = numel(xi);

if ~exist('bandwidth', 'var') || isempty(bandwidth)
  sig = median(abs(xi - median(xi)))/0.6745;
  if sig <= 0, sig = max(xi) - min(xi); end
  if sig > 0
    bandwidth = sig * (1/nframe)^(1/6);
  else
    bandwidth = 1;
  end

  fprintf('bandwidth: %f\n', bandwidth);
end

if ~exist('box', 'var') || isempty(box)
  is_box = false;
else
  is_box = true;
end

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nframe, 1);
  weight = weight./sum(weight);
end

%% compute the kernel density estimates
if is_box
  dx = bsxfun(@minus, data(:, 1), xi);
  dx = dx - round(dx./box(1))*box(1);
  dx2 = (dx./bandwidth(1)).^2;
else
  dx2 = (bsxfun(@minus, data(:, 1), xi)./bandwidth(1)).^2;
end
dx2 = exp(-0.5 * dx2)./(sqrt(2*pi).*bandwidth(1));
f = sum(bsxfun(@times, dx2, weight));
%f = logsumexp2(-0.5*(bsxfun(@times, dx2, log(weight))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = logsumexp2(x)
max_x = max(x);
exp_x= exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x)) + max_x;

