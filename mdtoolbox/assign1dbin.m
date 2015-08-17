function [index, histogram] = assign1dbin(data, edge)
%% assign1dbin
% assign 1-D bin index to the given 1-D data by binning their values
%
%% Syntax
%# index = assign1dbin(data, edge);
%# [index, histogram] = assign1dbin(data, edge)
%#
%% Description
% This routine assigns bin index to the given data by binning
% their values. 
%
% For each data(i), bin index k is assgined if data(i) is
%   edge(k) ≤ data(i) < edge(k+1)
% 
% NaN is assigned to data outside the given edges.
% Also, data outside the given edges are not counted in the
% histogram. 
%
% * data - data to be assigned 
%          [double nframe x 1]
% * edge - edges of bins
%          [double nedge x 1]
%
% * index     - bined index data
%               [integer nframe x 1]
% * histogram - binned histogram
%               [integer nedge-1 x 1]
%
%% Example
%# data = rand(100000, 2);
%# index = assign1dbins(data(:, 1), 0:0.1:1);
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
% 
%% See also
% assign2dbin, assignvoronoi
%

%% setup
if isrow(data)
  data = data';
end
nedge = numel(edge);

%% construct histogram and binning
[histogram, index] = histc(data, edge);

%% eliminate samples outside the specified bin edges
index(index == nedge) = NaN;
index(index == 0) = NaN;
histogram(nedge) = [];

