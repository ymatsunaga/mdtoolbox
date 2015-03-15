function [pmf, xi, yi] = calcpmf2d(data, xi, yi, bandwidth, box, weight)
%% calcpmf2d
% calculate 2D potential of mean force from 2D-data by using kernel density estimator
%
%% Syntax
% [pmf, xi, yi] = calcpmf2d(data)
% pmf = calcpmf2d(data, xi, yi)
% pmf = calcpmf2d(data, xi, yi, bandwidth, box, weight)
%
%% Description
%
% * data      - 2-dimensional data
%               [nstep x 2 double array]
% * xi        - equally spaced grid in x-axis
%               [n x 1 or 1 x n double]
% * yi        - equally spaced grid in y-axis 
%               [m x 1 or 1 x m double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 2 double]
% * box       - PBC box size if data is periodic.
%               If not given or empty, data is considered to be non-periodic. 
%               [1 x 2 double]
% * weight    - weights of observables in data. 
%               By default, uniform weight is used. 
%               [nstep x 1 double]
% * pmf       - 2-dimensional potential of mean force
%               [m x n double array]
%
%% Example
%# xi = -180:2:180; % grids in x-axis
%# yi = -180:2:180; % grids in y-axis
%# pmf = calcpmf2d([phi psi], xi, yi, [3.0 3.0], [360 360]);
%# s   = getconstants(); % get Boltzmann constant in kcal/mol/K
%# T   = 300.0;          % set temperature
%# pmf = s.KB*T*pmf;     % convert unit from KBT to kcal/mol
%# landscape(xi, yi, pmf, 0:0.25:6);
%
%% See also
% calcpmf, ksdensity2d
%

%% setup
if ~exist('xi', 'var')
  xi = [];
end

if ~exist('yi', 'var')
  yi = [];
end

if ~exist('bandwidth', 'var')
  bandwidth = [];
end

if ~exist('box', 'var')
  box = [];
end

if ~exist('weight', 'var')
  weight = [];
end

%% compute 2-dimensional potential of mean force
pdf = ksdensity2d(data, xi, yi, bandwidth, box, weight);
pdf(pdf < realmin('single')) = NaN;
pmf = -log(pdf);
pmf = pmf - min(pmf(:));
pmf = pmf';

