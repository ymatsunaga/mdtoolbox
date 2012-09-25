function [bin_kn, data0] = assign1dbins(data_kn, nbins, data_subsampled)
%% assign1dbins
%
%% Syntax
%#
%
%% Description
%
%% Example
%#
% 
%% See also
%
%% References
%

data_min = min(data_subsampled);
data_max = max(data_subsampled);
%data_min = 28;
%data_max = 40;
edges = linspace(data_min, data_max, nbins+1);
edges(end) = inf;
[h, bin_kn] = histc(data_kn, edges);
edges(end) = [];
data0 = edges + 0.5*(edges(2) - edges(1));

