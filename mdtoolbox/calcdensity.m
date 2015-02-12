function density = calcdensity(trj, xi, yi, zi, bandwidth, weight)
%% calcdensity
% compute 3-dimensional density of atoms from trajectory data
%
%% Syntax
%# density = calcdensity(trj);
%# density = calcdensity(trj, edges);
%
%% Description
%
% * trj       - trajectory data [nstep x natom3 double]
% * xi        - equally spaced grid in the x-axis on which
%               the density is estimated [1 x nx double]
% * yi        - equally spaced grid in the y-axis on which
%               the density is estimated [1 x ny double]
% * zi        - equally spaced grid in the z-axis on which
%               the density is estimated [1 x ny double]
% * bandwidth - bandwidth of the Gaussian kernel function
%               [1 x 3 double]
% * weight    - weights of trajectory snapshots. 
%               By default, uniform weight is used. 
%               [nstep x 1 double]
% * density   - volumetric data of 3-dimensional density
%               [nx x ny x nz double]
%               
%% Example:
% [trj, box] = readdcd('run.dcd');
% xi = linspace(-10, 10, 100);
% yi = linspace(-10, 10, 100);
% zi = linspace(-10, 10, 100);
% density = density(trj, xi, yi, zi, [1.0 1.0 1.0]);
% writedx('density.dx', density, xi, yi, zi);
%

%% preparation
nstep  = size(trj, 1);
natom3 = size(trj, 2);
natom  = natom3./3;
%bulk_density = (natom./(nbins.^3));

nx = numel(xi);
ny = numel(yi);
nz = numel(zi);
density = zeros(nx, ny, nz);

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nstep, 1);
end
weight = weight./sum(weight);

%% compute density
for istep = 1:nstep
  data = reshape(trj(istep, :), 3, natom)';
  f = ksdensity3d(data, xi, yi, zi, bandwidth);
  density = density + f*natom*weight(istep);
end
%density = density./bulk_density;
%density = density./nstep;

