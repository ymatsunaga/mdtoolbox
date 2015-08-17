function [index, center_x, center_y, edge_x, edge_y] = assign2dbins(data, edge_x, edge_y)
%% assign2dbins
% assign 1-dimensional indices to the 2-dimensional data
%
%% Syntax
%# index = assign2dbins(data);
%# index = assign2dbins(data, edge_x, edge_y);
%# [index, center_x, center_y] = assign2dbins(data, edge_x, edge_y);
%#
%# index = assign2dbins(data, nbin_x, nbin_y);
%# [index, center_x, center_y, edge_x, edge_y] = assign1dbins(data, nbin)
%
%% Description
% This routine assigns 1-dimensional bin indices to the given data
% in 2-dimensional space. 
%
% * data     - data to be assigned 
%              [double nframe x 1]
% * edge_x   - edges for bins in the 1st dimension
%              [double nframe+1 x 1]
% * edge_y   - edges for bins in the 2nd dimension
%              [double nframe+1 x 1]
% * index    - indices of bins.
%              [integer nframe x 1]
% * center_x - centers of bins in the 1st dimension
%              [double nframe x 1]
% * center_y - centers of bins in the 2nd dimension
%              [double nframe x 1]
% * nbin_x   - the number of bins in the 1st dimension
%              [integer scalar]
% * nbin_y   - the number of bins in the 2nd dimension
%              [integer scalar]
%
%% Example
%# data = rand(100000, 2);
%# [index, center, edge] = assign2dbins(data, 0:0.2:1, 0:0.5:1);
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
% 
%% See also
% assign1dbins, assignvoronoi
%

if isrow(data)
  data = data';
end

if (~exist('edge_x', 'var') || isempty(edge_x))
  edge_x = 10;
end

if (~exist('edge_y', 'var') || isempty(edge_y))
  edge_y = 10;
end

if numel(edge_x) == 1
  nbin_x = edge_x;
  data_min = min(data(:, 1));
  data_max = max(data(:, 1));
  edge_x = linspace(data_min, data_max + nbin_x*eps, nbin_x+1);
else
  nbin_x = numel(edge_x) - 1;
end

if numel(edge_y) == 1
  nbin_y = edge_y;
  data_min = min(data(:, 2));
  data_max = max(data(:, 2));
  edge_y = linspace(data_min, data_max + nbin_y*eps, nbin_y+1);
else
  nbin_y = numel(edge_y) - 1;
end

[index_x, center_x] = assign1dbins(data(:, 1), edge_x);
[index_y, center_y] = assign1dbins(data(:, 2), edge_y);
index = nbin_y*(index_x-1) + index_y;

