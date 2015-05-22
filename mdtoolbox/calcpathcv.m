function [progress, distance] = calcpathcv(path, data, lambda)
%% calcpathcv
% calculate the progress of the pathway as well as the distance from the closest point along the pathway
%
%% Syntax
%# [progress, distance] = calcpathcv(path, data);
%# [progress, distance] = calcpathcv(path, data, lambda);
%
%% Description
% This routine calculates the progress of the pathway as well as
% the distance from the point along the pathway. 
%
% * path     - path way data consists of m-points defined in n-dimensional space
%              [m x ndim double]
% * data     - some trajectory or data set to be evaluated
%              [nstep x ndim double]
% * progress - progress of the pathway
%              [nstep x 1 double]
% * distance - distance from the closest point along the pathway
%              [nstep x 1 double]
%
%% Example
%# data = rand(100000, 2);
%# [progress, distance] = calcpathcv([0.3 0.5; 0.6 0.5; 0.9 0.5], data);
%# scatter(data(:, 1), data(:, 2), 20, progress, 'filled');
%# scatter(data(:, 1), data(:, 2), 20, distance, 'filled');
% 
%% References
% D. Branduardi, F. L. Gervasio, and M. Parrinello, 
% J. Chem. Phys. 126, 054103 (2007).
%
%% See also
% assignvoronoi
%

%% setup
[m, ndim1] = size(path);
[nstep, ndim2] = size(data);
assert(ndim1 == ndim2, 'dimensions of path and data do not match...');
ndim = ndim1;

%% calculation
progress = zeros(nstep, 1);
distance = zeros(nstep, 1);
d = zeros(nstep, m);

%% determine lambda parameter
d_pathpoints = 0;
for i = 1:(m-1)
  d_pathpoints = d_pathpoints + sum((path(i,:) - path(m-1,:)).^2);
end
d_pathpoints = d_pathpoints/(m-1);
lambda = 2.3/d_pathpoints;
%lambda = 1.5/d_pathpoints;

%% calculate CVs
for i = 1:m
  dev = bsxfun(@minus, path(i, :), data);
  m
  whos dev
  d(:, i) = sum(dev.^2, 2);
end

if nargout > 0
  for istep = 1:nstep
    progress(istep, 1) = sum((1:m).*exp(-lambda*d(istep, :)), 2) ./ sum(exp(-lambda*d(istep, :)), 2);
  end
end

if nargout > 1
  for istep = 1:nstep
    distance(istep, 1) = - (1/lambda) * logsumexp(-lambda*d(istep, :));
  end
end

%% logsumexp (input should be a vector)
function s = logsumexp(x)
max_x = max(x);
exp_x = exp(x - max_x);
s = log(sum(exp_x)) + max_x;

