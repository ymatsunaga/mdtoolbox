function [z, xi, yi] = calchistpmf2d(x, nbin)
% calculate 2D potential of mean force from the scattered 2D-data (using histogram)
%
% function [z,xi,yi] = fene2d_hist(x,nbin)
%
% %AKの粗視化モデルの主成分データを読み込む
% load p.mat;
% %ヒストグラム計算に用いるbinの数
% nbin = 100;
% [z,xi,yi] = fene2d_hist(p(:,[1 2]),nbin);
% imagesc(xi,yi,z); axis xy;
% xlabel('PC1','FontSize',45); ylabel('PC2','FontSize',45); colorbar; plot_format
%

xd = x(:,1);
yd = x(:,2);

if nargin == 1
  nbin = 100;
end

[z,xi,yi] = hist2d([xd yd],nbin);
%z = z ./ bin_width.^2;

z = -log(z);
z_max = max(max(-z));
z = z + z_max;

