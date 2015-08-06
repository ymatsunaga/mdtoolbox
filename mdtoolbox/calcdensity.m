function [density, xi, yi, zi] = calcdensity(trj, xi, yi, zi, bandwidth, box, weight)
%% calcdensity
% compute 3-dimensional density of atoms from trajectory data
%
%% Syntax
%# [density, xi, yi, zi] = calcdensity(trj);
%# density = calcdensity(trj, xi, yi, zi);
%# density = calcdensity(trj, xi, yi, zi, bandwidth);
%# density = calcdensity(trj, xi, yi, zi, bandwidth, weight);
%
%% Description
%
% * trj       - trajectory data [nframe x natom3 double]
% * xi        - uniformly-spaced grid in the x-axis 
%               on which the density is estimated [1 x nx double]
% * yi        - uniformly-spaced grid in the y-axis 
%               on which the density is estimated [1 x ny double]
% * zi        - uniformly-spaced grid in the z-axis 
%               on which the density is estimated [1 x ny double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 3 double]
% * weight    - weights of trajectory snapshots. 
%               By default, uniform weight is used. 
%               [nframe x 1 double]
% * density   - volumetric data of 3-dimensional density
%               Note that the density value is multiplied by the volumetric element.
%               Thus, the value of this matrix is actually not
%               density but the average number of atoms.
%               [nx x ny x nz double]
%               
%% Example:
% [trj, box] = readdcd('run.dcd');
% psf = readpsf('run.psf');
% index = selectname(psf.residue_name, 'TIP3');
% [density, xi, yi, zi] = calcdensity(trj(:, to3(index)));
% isosurface(xi, yi, zi, density, 0.5);
% writedx('density.dx', density, xi, yi, zi);
% $ vmd -dx density.dx
%

%% preparation
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom  = natom3./3;
%bulk_density = (natom./(nbins.^3));

if ~exist('xi', 'var') || isempty(xi)
  g_min = min(trj(1, 1:3:end));
  g_max = max(trj(1, 1:3:end));
  xi = g_min:1:g_max;
  disp(sprintf('message: xi = %f:1.0:%f is used.', g_min, g_max));
end

if ~exist('yi', 'var') || isempty(xi)
  g_min = min(trj(1, 2:3:end));
  g_max = max(trj(1, 2:3:end));
  yi = g_min:1:g_max;
  disp(sprintf('message: yi = %f:1.0:%f is used.', g_min, g_max));
end

if ~exist('zi', 'var') || isempty(zi)
  g_min = min(trj(1, 3:3:end));
  g_max = max(trj(1, 3:3:end));
  zi = g_min:1:g_max;
  disp(sprintf('message: zi = %f:1.0:%f is used.', g_min, g_max));
end

if ~exist('bandwidth', 'var') || isempty(bandwidth)
  bandwidth = [1.0 1.0 1.0];
  disp('message: bandwidth = [1.0 1.0 1.0] is used.');
end

if ~exist('box', 'var')
  box = [];
else (size(box, 1) == 1)
  box = repmat(box, nframe, 1);
end

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nframe, 1);
  disp('message: weights is considered as uniform.');
end
weight = weight./sum(weight);

%% compute density
nx = numel(xi);
ny = numel(yi);
nz = numel(zi);
density = zeros(nx, ny, nz);
f = zeros(nx, ny, nz);

if isempty(box)
  for iframe = 1:nframe
    iframe
    data = reshape(trj(iframe, :), 3, natom)';
    f = ksdensity3d(data, xi, yi, zi, bandwidth, []);
    density = density + f*weight(iframe);
  end
else
  for iframe = 1:nframe
    iframe
    data = reshape(trj(iframe, :), 3, natom)';
    f = ksdensity3d(data, xi, yi, zi, bandwidth, box(iframe, :));
    density = density + f*weight(iframe);
  end
end
density = density*(natom./(abs(xi(2)-xi(1))*abs(yi(2)-yi(1))*abs(zi(2)-zi(1))));
%density = density./bulk_density;

