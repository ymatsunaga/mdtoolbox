function [progress, distance] = calcpathcv(path, data)
%% calcpathcv
% calculate the progress of the pathway and the distance from the closest point along the pathway
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
%              [nframe x ndim double]
% * progress - progress of the pathway
%              [nframe x 1 double]
% * distance - distance from the closest point along the pathway
%              [nframe x 1 double]
%
%% Example
%# path = [0.25 0.3; 0.5 0.6; 0.75 0.3];
%# data = rand(100000, 2);
%# [progress, distance] = calcpathcv(path, data);
%# subplot(1, 2, 1);
%# scatter(data(:, 1), data(:, 2), 20, progress, 'filled');
%# hold on;
%# scatter(path(:, 1), path(:, 2), 200, 'k', 'filled');
%# axis square; 
%# subplot(1, 2, 2);
%# scatter(data(:, 1), data(:, 2), 20, distance, 'filled');
%# hold on;
%# scatter(path(:, 1), path(:, 2), 200, 'k', 'filled');
%# axis square;
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
[nframe, ndim2] = size(data);
assert(ndim1 == ndim2, 'dimensions of path and data do not match...');

%% calculation
d = zeros(nframe, m);

% determine lambda parameter
d_pathpoints = 0;
for i = 1:(m-1)
  d_pathpoints = d_pathpoints + sqrt(sum((path(i+1,:) - path(i,:)).^2));
end
d_pathpoints = d_pathpoints/(m-1);
lambda = 2.3/(d_pathpoints.^2);
%lambda = 1.5/d_pathpoints;

% calculate pathway CVs
for i = 1:m
  dev = bsxfun(@minus, path(i, :), data);
  d(:, i) = sum(dev.^2, 2);
end

log_m = log(1:m);
log_numerator   = logsumexp_array(bsxfun(@minus, log_m, lambda*d));
log_denominator = logsumexp_array(-lambda*d);
progress = exp(log_numerator - log_denominator);

distance = - (1/lambda) * logsumexp_array(-lambda*d);

%% logsumexp_array (input should be an array matrix, summed over the 2nd dimension)
function s = logsumexp_array(x)
max_x = max(x, [], 2);
exp_x = exp(bsxfun(@minus, x, max_x));
s = log(sum(exp_x, 2)) + max_x;

