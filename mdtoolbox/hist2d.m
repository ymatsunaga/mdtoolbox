function [z,xi,yi] = hist2d(x,nbin)
% calculate free-energy surface
%
% function [z,xi,yi] = hist2d(x,nbin)
%
% %AKの粗視化モデルの主成分データを読み込む
% load p.mat;
% %ヒストグラム計算に用いるbinの数
% nbin = 100;
% [z,xi,yi] = hist2d(p(:,[1 2]),nbin);
% imagesc(xi,yi,z); axis xy;
% xlabel('PC1','FontSize',45); ylabel('PC2','FontSize',45); colorbar; plot_format
%

xd = x(:,1);
yd = x(:,2);

if nargin == 1
  nbin = 100;
end  

% calc 2-d histogram of x(:,1:2)
% see http://blogs.mathworks.com/videos/2010/01/22/advanced-making-a-2d-or-3d-histogram-to-visualize-data-density/
n  = nbin;
xi = linspace(min(xd(:)),max(xd(:)),n);
yi = linspace(min(yd(:)),max(yd(:)),n);
xr = interp1(xi,1:numel(xi),xd,'nearest')';
yr = interp1(yi,1:numel(yi),yd,'nearest')';
z = accumarray([xr' yr'], 1, [n n]);
z = z';

