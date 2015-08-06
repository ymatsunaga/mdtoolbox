function [index, dmin] = assignvoronoi(ref, data)
%% assignvoronoi
% assign Voronoi cell indices to the given data by checking the distances between ref and data
%
%% Syntax
%# index = assignvoronoi(ref, data)
%# [index, dmin] = assignvoronoi(ref, data)
%
%% Description
% This routine assigns Voronoi cell indices to the elements of data
% by checking the distances between ref and data. For the distance,
% n-dimensional Euclid distance is used. 
% Data points of index=i belong to the Voronoi cell of the
% reference point (or generating point) point ref(i,:).
%
% * ref   - generating points of Voronoi cells 
%           [n x ndim double]
% * data  - some trajectory or data set to be assigned 
%           [nframe x ndim double]
% * index - indices of Voronoi cells. indices correspond to the rows of ref 
%           [nframe x 1 double]
% * dmin  - distance between the data to the nearest generating point 
%           [nframe x 1 double]
%
%% Example
%# ref  = rand(10, 2);
%# data = rand(10000, 2);
%# index = assignvoronoi(ref, data);
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
%# hold on
%# scatter(ref(:, 1), ref(:, 2), 100, 'ro', 'filled');
% 
%% See also
% assign1dbins, assign2dbins, 
%

%% setup
[~, ndim1] = size(ref);
[nframe, ndim2] = size(data);
assert(ndim1 == ndim2, 'dimensions of ref and data do not match...');

%% calculation
index = zeros(nframe, 1);
dmin = zeros(nframe, 1);

for iframe = 1:nframe
  dev = bsxfun(@minus, ref, data(iframe, :));
  dist = sum(dev.^2, 2);
  [dmin_each, index_each] = min(dist);
  index(iframe) = index_each;
  dmin(iframe) = sqrt(dmin_each);
end

