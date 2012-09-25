function [bin_kn_2d,data0,nbins_found,bin_indices] = assign2dbins(data_kn,nbins,data_subsampled,data_subsampled_2d)
%% assign2dbins
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

[K, N_max] = size(data_kn);
data_min = min(data_subsampled);
data_max = max(data_subsampled);
%data_min = 28;
%data_max = 40;
edges = linspace(data_min, data_max, nbins+1);
edges(end) = inf;
[h, bin_kn] = histc(data_kn, edges);

bin_kn_2d = bin_kn;
nbins_found = 0;
bin_indices = [];
for k = 1:K
  for n = 1:nbins
    %in_bin = (bin_kn(k,:) == n);
    [h, bin_subsampled] = histc(data_subsampled_2d{k}, edges);
    in_bin = (bin_subsampled == n);
    bin_count = sum(in_bin);
    if (bin_count > 0)
      nbins_found = nbins_found + 1;
      in_bin = (bin_kn(k, :) == n);
      bin_kn_2d(k, in_bin) = nbins_found;
      bin_indices = [bin_indices; [k n]];
    end
  end
end

edges(end) = [];
data0 = edges + 0.5*(edges(2) - edges(1));

