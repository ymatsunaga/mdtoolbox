function [index, histogram] = assign2dbin(data, edge_x, edge_y)
%% assign2dbin
% assign 2-D bin index to the given 2-D data by binning their values
%
%% Syntax
%# index = assign2dbin(data, edge_x, edge_y);
%# [index, histogram] = assign2dbin(data, edge_x, edge_y)
%
%% Description
% This routine assigns bin index to the given data by binning
% their values. 
%
% For each data(i, :), bin index k1 and k2 is assgined if data(i, :) is
%   edge_x(k1) ≤ data(i, 1) < edge_x(k1+1)
%   edge_y(k2) ≤ data(i, 2) < edge_y(k2+1)
% 
% NaN is assigned to data outside the given edges.
% Also, data outside the given edges are not counted in the
% histogram. 
%
% * data     - 2-D data to be assigned 
%              [double nframe x 2]
% * edge_x   - edges for bins in the 1st dimension
%              [double (nedge_x) x 1]
% * edge_y   - edges for bins in the 2nd dimension
%              [double (nedge_y) x 1]
%
% * index     - bined 2-D index data
%               [integer nframe x 2]
% * histogram - binned histogram
%               [integer (nedge_x-1) x (nedge_y-1)]
%
%% Example
%# data = rand(100000, 2);
%# index = assign2dbin(data, 0:0.2:1, 0:0.5:1);
%# scatter(data(:, 1), data(:, 2), 5, index(:, 1) + 5*(index(:, 2)-1), 'filled');
% 
%% See also
% assign1dbin, assignvoronoi
%

%% setup
nedge_x = numel(edge_x);
nedge_y = numel(edge_y);

%% assing bin-index to samples
[~, index_x] = histc(data(:, 1), edge_x);
[~, index_y] = histc(data(:, 2), edge_y);

%% eliminate samples outside the specified bin edges
id1 = (index_x == nedge_x);
id2 = (index_x == 0);
id12 = id1 | id2;
index_x(id12) = NaN;

id3 = (index_y == nedge_y);
id4 = (index_y == 0);
id34 = id3 | id4;
index_y(id34) = NaN;

index = [index_x index_y];

%% construct histogram with accumulation
if nargout > 1
  id = id12 | id34;
  index = index(~id, :);
  z = accumarray(index, 1, [(nedge_x-1) (nedge_y-1)]);
end

