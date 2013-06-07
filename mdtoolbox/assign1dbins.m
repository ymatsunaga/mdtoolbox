function [index, center, edge] = assign1dbins(data, nbin)
%% assign1dbins
% assign bin indices to the given data by binning their values
%
%% Syntax
%# index = assign1dbins(data, nbin);
%# [index, center, edge] = assign1dbins(nbin, data)
%#
%# index = assign1dbins(data, edge);
%# [index, center] = assign1dbins(nbin, edge)
%
%% Description
% This routine assigns bin indices to the given data by binning
% their values. Bins are numbered in acending order of the data
% values. For example min(data) belongs to the 1st bin, and
% max(data) belongs to the last bin.
%
% * nbin   - the number of bins
%            [integer scalar]
% * data   - some trajectory or data set to be assigned 
%            [double nstep x 1]
% * index  - indices of bins.
%            [integer nstep x 1]
% * center - centers of bins
%            [double nstep x 1]
% * edge   - edges of bins
%            [double nstep+1 x 1]
%
%% Example
%# data = rand(10000, 2);
%# [index, center, edge] = assign1dbins(10, data(:, 1));
%# scatter(data(:, 1), data(:, 2), 5, index, 'filled');
% 
%% See also
% assign2dbins, assignvoronoi
%

if isrow(data)
  data = data';
end

if (nargin < 2) || (numel(nbin) == 0)
  nbin = 100;
end

if numel(nbin) == 1
  data_min = min(data(:, 1));
  data_max = max(data(:, 1));
  edge = linspace(data_min, data_max + nbin*eps, nbin+1);
else
  edge = nbin;
end
clear nbin;

[~, index] = histc(data, edge);
center = edge + 0.5*(edge(2) - edge(1));
center(end) = [];

