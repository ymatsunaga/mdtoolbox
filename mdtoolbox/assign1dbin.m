function [index, center, edge] = assign1dbins(data, edge)
%% assign1dbins
% assign bin indices to the given data by binning their values
%
%% Syntax
%# index = assign1dbins(data, edge);
%# [index, center] = assign1dbins(data, edge)
%#
%# index = assign1dbins(data, nbin);
%# [index, center, edge] = assign1dbins(data, nbin)
%
%% Description
% This routine assigns bin indices to the given data by binning
% their values. Bins are numbered in acending order of the data
% values. For example min(data) belongs to the 1st bin, and
% max(data) belongs to the last bin.
%
% * data   - data to be assigned 
%            [double nframe x 1]
% * edge   - edges of bins
%            [double nframe+1 x 1]
% * index  - indices of bins.
%            [integer nframe x 1]
% * center - centers of bins
%            [double nframe x 1]
% * nbin   - the number of bins
%            [integer scalar]
%
%% Example
%# data = rand(100000, 2);
%# [index, center] = assign1dbins(data(:, 1), 0:0.1:1);
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
% 
%% See also
% assign2dbins, assignvoronoi
%

if isrow(data)
  data = data';
end

if ~exist('edge', 'var') || isempty(edge)
  edge = 10;
end

if numel(edge) == 1
  nbin = edge;
  data_min = min(data(:, 1));
  data_max = max(data(:, 1));
  edge = linspace(data_min, data_max + nbin*eps, nbin+1);
else
  nbin = numel(edge) - 1;
end

[~, index] = histc(data, edge);
center = 0.5*(edge(2:end) + edge(1:(end-1)));

