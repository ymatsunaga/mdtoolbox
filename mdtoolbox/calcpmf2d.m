function [z,xi,yi] = calcpmf2d(x,nbin)
%% calcpmf2d
% calculate 2D potential of mean force from the scattered 2D-data (using kernel density estimator)
%
% function [z,xi,yi] = fene2d(x,nbin)
%
% %AKの粗視化モデルの主成分データを読み込む
% load p.mat;
% %ヒストグラム計算に用いるbinの数
% nbin = 100;
% [z,xi,yi] = fene2d(p(:,[1 2]),nbin);
% imagesc(xi,yi,z); axis xy;
% xlabel('PC1','FontSize',45); ylabel('PC2','FontSize',45); colorbar; plot_format
%

xd = x(:,1);
yd = x(:,2);

if ~exist('nbin', 'var') || isempty(nbin)
  nbin = 100;
end

xi = linspace(min(xd),max(xd),nbin);
yi = linspace(min(yd),max(yd),nbin);
[bandwidth,z,xi,yi]=kde2d([xd yd],nbin);
xi = xi(1,:);
yi = yi(:,1)';

z(z<0) = 0.0;
z = -log(z);
%z = -log(abs(z));
z_max = max(max(-z));
z = z + z_max;

