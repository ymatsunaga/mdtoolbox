function gobj = plotcolor(x, y, c, varargin)
%% plotcolor
% plot colored lines according to the input vector c
%
%% Syntax
%# plotcolor(x, y, c)
%# gobj = plotcolor(x, y, c, LineseriesProperties)
%
%% Description
%
% * x        - values in the x-axis [n x 1 double]
% * y        - values in the y-axis [n x 1 double]
% * c        - values for the line color [n x 1 double]
% * varargin - Lineseries Properties such as 'LineWidth', 2.
% * gobj     - graphics object (formerly called as 'handle graphics')
%
%% Example
%# parm = readparm('ak.parm');
%# trj = readnetcdf('ak.nc');
%# index_ca = find(selectname(parm.atom_name, 'CA'));
%# indexOfCluster = clusteringbykmeans(trj(:, to3(index_ca)), 3);
%# rmsd = superimpose(trj(1, :), trj, index_ca);
%# plotcolor(1:numel(rmsd), rmsd, indexOfCluster, 'linewidth', 2);
%

%% preparation
nframe1 = numel(x);
nframe2 = numel(y);
nframe3 = numel(c);
assert(nframe1 == nframe2, 'x and y must be the same size');
assert(nframe1 == nframe3, 'x and c must be the same size');
nframe = nframe1;

%% determine RGB values from c vector
% transform c to [0, 1]
c_min = min(c);
c_max = max(c);
c = c - c_min;
c = c./max(c);

% transform [0, 1] to [1 size(cmap, 1)]
cmap = colormap;
ncolor = size(cmap, 1);
c = c * (ncolor-1) + 1;

% interpolate RGB values by using cubic spline 
cs = spline(1:ncolor, cmap');
rgb = ppval(cs, c);
rgb = rgb';

rgb(rgb < 0) = 0;
rgb(rgb > 1) = 1;

%% plot lines
if numel(varargin) < 2
  for i = 2:nframe
    gobj = line([x(i-1) x(i)], [y(i-1) y(i)], 'Color', rgb(i-1, :));
  end
else
  for i = 2:nframe
    gobj = line([x(i-1) x(i)], [y(i-1) y(i)], 'Color', rgb(i-1, :), varargin{:});
  end
end

caxis([c_min c_max]);





