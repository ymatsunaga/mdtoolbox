function rc = exportas(basename)
%% exportas
% export fig, eps, png, tiff files of the current figure
%
%% Syntax
% exportas(basename)
%
%% Description
% create fig, eps, png, tiff files with this single function.
% For example, exportas('plot') creates the following files:
% plot.fig
% plot.eps
% plot.png
% plot.tiff
% 
% * basename  - basename of output files [chars]
%         
%% Example
% plot(randn(1000,1))
% formatplot;
% xlabel('step','FontSize',50)
% ylabel('value','FontSize',50)
% exportas('myplot');
%
%% See also
% formatplot
% 

f = [basename '.fig'];
saveas(gcf,f);

set(gcf,'PaperPositionMode','auto');

f = [basename '.eps'];
print('-depsc',f);

f = [basename '.png'];
print('-dpng',f);

f = [basename '.tiff'];
print('-dtiff',f);

