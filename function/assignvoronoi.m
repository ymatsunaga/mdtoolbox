function [index, dmin] = assignvoronoi(ref, data)
%% assignvoronoi
% assign Voronoi cell indices data by checking the distances between ref and data
%
%% Syntax
%# index = assignvoronoi(ref, data)
%# [index, dmin] = assignvoronoi(ref, data)
%
%% Description
% This routine assigns Voronoi cell indices to the elements of data
% by checking the distances between ref and data. For the distance,
% n-dimensional Euclid distance is used. 
%
% * ref   - generating points of Voronoi cells [n x ndim double]
% * data  - some trajectory or data set to be assigned [nstep x ndim double]
% * index - indices of Voronoi cell. indices correspond to the rows of ref [nstep double]
% * box   - distance between the data to the nearest generating point [nstep double]
%
%% Example
%# ref  = rand(10, 2);
%# data = rand(10000, 2);
%# [index, dmin] = assignvoronoi(ref, data);
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
%# hold on
%# scatter(ref(:, 1), ref(:, 2), 100, 'ro', 'filled');
% 
%% See also
% assign1dbins, assign2dbins, 
%

%% setup
[ncell, ndim1] = size(ref);
[nstep, ndim2] = size(data);
assert(ndim1 == ndim2, 'dimensions of ref and data do not match...');

%% calculation
index = zeros(nstep, 1);
dmin = zeros(nstep, 1);

for istep = 1:nstep
  dev = bsxfun(@minus, ref, data(istep, :));
  dist = sum(dev.^2, 2);
  [dmin_each, index_each] = min(dist);
  index(istep) = index_each;
  dmin(istep) = sqrt(dmin_each);
end


