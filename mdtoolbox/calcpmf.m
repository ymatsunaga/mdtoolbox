function [z,xi] = calcpmf(x,nbin)
%% calcpmf
% calculate 1D potential of mean force from the scattered 1D-data (using kernel density estimator)
%
% function [z,xi] = fene1d(x,nbin)
%
% %AKの粗視化モデルの主成分データを読み込む
% load p.mat
% %ヒストグラム計算に用いるbinの数
% nbin = 1000;
% [z,xi] = fene1d(p(:,1),nbin);
% plot(xi,z,'-'); xlabel('PC1','FontSize',45); ylabel('\Delta {\itF} ({\itk_{B}T})','FontSize',45); plot_format;
%

x = x(:,1);
nstep = size(x,1);

if nargin == 1
  nbin = 100;
end

xi = linspace(min(x),max(x),nbin);
z = ksdensity(x,xi);

z = -log(z);
z_max = max(max(-z));
z = z + z_max;

