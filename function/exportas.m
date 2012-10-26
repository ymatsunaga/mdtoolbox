function rc = exportas(basename, style)
%% exportas
% export fig, eps, png, tiff files of the current figure
%
%% Syntax
% exportas(basename)
% exportas(basename, style)
%
%% Description
% This code creates fig, eps, png, compressed tiff images.
% For example, exportas('plot') creates the following images:
% plot.fig
% plot.eps
% plot.png
% plot.tiff
%
% By default, these image are created in WYSWIG style, 
% i.e, you get the same pixels and sizes as the figure window.
% If style option (either 'single' or 'double') is given
% the width of the images are adjusted for 
% single-column width (8.3cm) or double-column width (17.35cm), 
% and the resolution becomes 300 ppi for publication.
% 
% * basename  - basename of output files [chars]
% * style     - speficy the image style. 'single' or 'double' [chars]
%         
%% Example
% plot(randn(1000, 1))
% formatplot;
% xlabel('step', 'FontSize', 30)
% ylabel('value', 'FontSize', 30)
% exportas('myplot');
% exportas('myplot_singlecolumn', 'single');
% exportas('myplot_doublecolumn', 'double');
%
%% See also
% formatplot
% 
%% References
% Guidelines for Figure Preparation (PLoS Comput Biol)
% http://www.ploscompbiol.org/static/figureGuidelines.action
% 

%% setup
ppi = 300;

if nargin < 2
  style = 'wyswig';
end

if strncmpi(style, 'single', numel('single'))
  image_width = 8.3;
elseif strncmpi(style, 'double', numel('double'))
  image_width = 17.35;
else
  error(sprintf('sorry, does not support your style option %s', style));
end
  
%% create figures
f = [basename '.fig'];
saveas(gcf, f);

if strncmpi(style, 'default', numel('single'))
  % WYSWIG
  set(gcf, 'PaperPositionMode', 'auto');
  
  f = [basename '.eps'];
  print('-depsc', f);
  
  f = [basename '.png'];
  print('-dpng', f);
  
  f = [basename '.tiff'];
  print('-dtiff', f);

else
  % single-column width with 300 ppi
  set(gcf, 'PaperPositionMode', 'manual');
  
  set(gcf, 'Units', 'centimeters');
  set(gcf, 'PaperUnits', 'centimeters');
  position      = get(gcf, 'Position');
  figure_width  = position(3);
  figure_height = position(4);
  image_height  = figure_height * (image_width./figure_width)
  
  set(gcf, 'Papersize', [image_width image_height]);
  set(gcf, 'PaperPosition', [0 0 image_width image_height]);
  
  f = [basename '.eps'];
  print('-depsc', sprintf('-r%d', ppi), f);
  
  f = [basename '.png'];
  print('-dpng', sprintf('-r%d', ppi), f);
  
  f = [basename '.tiff'];
  print('-dtiff', sprintf('-r%d', ppi), f);

end

