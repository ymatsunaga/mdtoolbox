function gobj = landscape(xi, yi, data, level)
%% landscape
% draws 2-dimensional free energy surface and its contour lines
%
%% Syntax
%# landscape(xi, yi, data)
%# gobj = landscape(xi, yi, data, level, level_max)
%
%% Description
%
% * xi    - regularly spaced grids in x-axis [n x 1 or 1 x n double vector]
% * yi    - regularly spaced grids in y-axis [m x 1 or 1 x m double vector]
% * data  - free energy data [m x n double array]
% * level - contour levels [doulbe vector]
% * gobj  - graphics object (formerly called as 'handle graphics')
%
%% Example
%#
%

%% preparation
data = data - min(data(:));
if ~exist('level', 'var') || isempty(level)
  level = 0:0.25:max(data(:));
end

level_max = max(level);

%% plot
data2 = data;
data2(data2 > level_max) = NaN;
gobj = pcolor(xi, yi, data2);
shading flat;
%colorbar;
axis([min(xi) max(xi) min(yi) max(yi)]);
axis xy;
formatplot2;

hold on;
gobj = contour(xi, yi, data2, level, 'linecolor', 'black');
hold off;

