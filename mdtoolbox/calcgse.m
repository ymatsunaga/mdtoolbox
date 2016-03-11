function [potential, density, energy, energy_total, force, xi, yi, zi] = calcgse(trj, charge, box, xi, yi, zi, sigma, weight, index_atom);
%% calcgse
% compute smooth (reciprocal) electrostatic energy, average potential/density using k-space Gaussian split Ewald
%
%% Syntax
%# potential = calcgse(trj, charge, box, xi, yi, zi);
%# potential = calcgse(trj, charge, box, xi, yi, zi, sigma);
%# potential = calcgse(trj, charge, box, xi, yi, zi, sigma, weight);
%# [potential, energy, energy_total, density] = calcgse(trj, charge, box, xi, yi, zi, sigma, weight);
%
%% Description
% This routine calculates (reciprocal part of) electrostatic energy.
% Specifically, electrostatic potential, energy, and smoothed
% charge density on the specified grid. These may be useful for visualization. 
%
% Input
% * trj    - trajectory [nframe x natom3 double]
% * charge - charge [natom x 1]
% * box    - PBC box size [nframe x 3 double]
% * xi     - equally spaced grid in the x-axis on which the potential/density is averaged [nx vector double]
% * yi     - equally spaced grid in the y-axis on which the potential/density is averaged [ny vector double]
% * zi     - equally spaced grid in the z-axis on which the potential/density is averaged [nz vector double]
% * sigma  - width of Gaussian function [scalar double]
% * weight - weights used for averaging. BY default, uniform weight is used [nframe vector double]
%
% Output
% * potential    - electrostatic potential on the specfied grid averaged over all frames [double, nx x ny nz]
% * density      - smoothed charge density on the specfied grid averaged over all frames [double, nx x ny nz]
% * energy       - electrostatic energy on the specfied grid averaged over all frames [double, nx x ny nz]
% * energy_total - time-series of total (reciprocal part of) electrostatic energy [double nframe]
%
%% Example
%# psf = readpsf('run.psf');
%# [trj, box] = readdcd('run.dcd');
%# xi = -32:32;
%# yi = -32:32;
%# zi = -32:32;
%# potential = calcgse(trj, psf.charge, box, xi, yi, zi);
%# writedx('potential.dx', potential, xi, yi, zi);
%
%% Reference
% Y. Shan, J. L. Klepeis, M. P. Eastwood, R. O. Dror, and D. E. Shaw, 
% J. Chem. Phys. 122, 54101 (2005).
% 

%% setup
coefficient = 332.0716;      % CHARMM
%coefficient = 332.05221729; % AMBER
%coefficient = 332.0636930;  % GROMACS

nframe  = size(trj, 1);
natom   = numel(charge);
energy  = [];
density = [];

if ~exist('trj', 'var') || isempty(trj)
  error('trj is required');
end

if ~exist('charge', 'var') || isempty(charge)
  error('charge is required');
end

if ~exist('box', 'var') || isempty(box)
  error('box is required');
end

if ~exist('xi', 'var') || isempty(xi)
  xi = linspace(-box(1, 1)/2, box(1, 1)/2, floor(box(1, 1))+1);
end

if ~exist('yi', 'var') || isempty(yi)
  yi = linspace(-box(1, 2)/2, box(1, 2)/2, floor(box(1, 2))+1);
end

if ~exist('zi', 'var') || isempty(zi)
  zi = linspace(-box(1, 3)/2, box(1, 3)/2, floor(box(1, 3))+1);
end

[xi_query, yi_query, zi_query] = meshgrid(xi, yi, zi);

if ~exist('sigma', 'var') || isempty(sigma)
  alpha = 0.34;
  sigma = 1./(alpha*sqrt(2));
  fprintf('message: sigma of %f A is used\n', sigma);
end

kappa = 0.35; % kappa should be in [0, 0.5]
%kappa = 0.45; % kappa should be in [0, 0.5]
%kappa = 0.50; % kappa should be in [0, 0.5]
sigma1 = sigma*sqrt(kappa);
sigma2 = sigma*sqrt(1 - 2*kappa);
bandwidth1 = [sigma1 sigma1 sigma1];

fprintf('message: sigma1 of %f A is used\n', sigma1);
fprintf('message: sigma2 of %f A is used\n', sigma2);

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nframe, 1);
  weight = weight./sum(weight);
end

if ~exist('index_atom') || isempty(index_atom)
  index_atom = 1:natom;
end

if islogical(index_atom)
  index_atom = find(index_atom);
end

if nargout >= 5
  force = zeros(natom, 3);
end

%% compute k-space Gaussian split Ewald
% pre-allocation
data = zeros(natom, 3);
potential = zeros(numel(yi), numel(xi), numel(zi));
density   = zeros(numel(yi), numel(xi), numel(zi));
energy    = zeros(numel(yi), numel(xi), numel(zi));
energy_total = zeros(nframe, 1);

for iframe = 1:nframe
  % determine grids for k-GSE
  xi_gse = linspace(-box(iframe, 1)/2, box(iframe, 1)/2, floor(box(iframe, 1))+1);
  yi_gse = linspace(-box(iframe, 2)/2, box(iframe, 2)/2, floor(box(iframe, 2))+1);
  zi_gse = linspace(-box(iframe, 3)/2, box(iframe, 3)/2, floor(box(iframe, 3))+1);
  xi_gse(end) = [];
  yi_gse(end) = [];
  zi_gse(end) = [];
  nx = numel(xi_gse);
  ny = numel(yi_gse);
  nz = numel(zi_gse);

  % 1: charge spreading by convolving Gaussian function in real-space
  data = reshape(trj(iframe, :), 3, natom)';
  d1 = ksdensity3d(data, xi_gse, yi_gse, zi_gse, bandwidth1, box(iframe, :), charge);
  d1 = permute(d1, [2,1,3]);
  d2 = d1;

  % FFT
  f = fft3d(d2);

  % 2: on-mesh charge spreading by convolution with Gaussian function exp(-sig^2*k^2/2)
  [X, Y, Z] = meshgrid(2*pi*ifftshift(-ceil((nx-1)/2):floor((nx-1)/2))./box(iframe, 1), ...
                     2*pi*ifftshift(-ceil((ny-1)/2):floor((ny-1)/2))./box(iframe, 2), ...
                     2*pi*ifftshift(-ceil((nz-1)/2):floor((nz-1)/2))./box(iframe, 3));
  gamma_mod = exp(-(sigma2.^2)*(X.^2 + Y.^2 + Z.^2)./2);

  % 3: solve poisson equation by convolution with Green function 4*pi/k^2
  gamma_mod = gamma_mod .* (4*pi) ./ (X.^2 + Y.^2 + Z.^2);
  gamma_mod(1, 1, 1) = 0;
  f = f.*gamma_mod;

  % 4: inverse FFT and get the on-mesh potential
  potential_gse = ifft3d(f);

  % 4': compute reciprocal energy
  d2 = d2.*potential_gse;
  d2 = coefficient * 0.5 * abs((xi_gse(2)-xi_gse(1)) * (yi_gse(2)-yi_gse(1)) * (zi_gse(2)-zi_gse(1))) * d2;
  energy_total(iframe) = sum(d2(:));

  % 4'': mesh interpolation to given fixed grids for averaging over frame
  %potential = potential + weight*interp3(xi_gse, yi_gse, zi_gse, potential_gse, xi_query, yi_query, zi_query, 'linear');
  potential = potential + weight(iframe)*interp3(xi_gse, yi_gse, zi_gse, potential_gse, xi_query, yi_query, zi_query, 'cubic');
  density   = density   + weight(iframe)*interp3(xi_gse, yi_gse, zi_gse, d1,            xi_query, yi_query, zi_query, 'cubic');
  energy    = energy    + weight(iframe)*interp3(xi_gse, yi_gse, zi_gse, d2,            xi_query, yi_query, zi_query, 'cubic');

  if nargout >= 5
    force_gse = calcforce(potential_gse, data, charge, sigma1, xi_gse, yi_gse, zi_gse, box(iframe, :), index_atom);
    force = force + weight(iframe)*force_gse;
  end
end

potential = permute(potential,[2,1,3]);
density   = permute(density,  [2,1,3]);
energy    = permute(energy,   [2,1,3]);

%%%%%% 3D FFT routines
function A_k = fft3d(A_r);
A_k = fft(fft(fft(A_r,[],1),[],2),[],3);

function A_r = ifft3d(A_k);
A_r = ifft(ifft(ifft(A_k,[],1),[],2),[],3,'symmetric');

%%%%%% force
function force = calcforce(potential, data, charge, sigma1, xi, yi, zi, box, index_atom);
nx = numel(xi);
ny = numel(yi);
nz = numel(zi);
natom = numel(charge);
natom_sub = numel(index_atom);

data = data(index_atom, :);

dx = bsxfun(@minus, data(:, 1), xi);
dx = dx - round(dx./box(1))*box(1);
dx2 = (dx./sigma1).^2;
dy = bsxfun(@minus, data(:, 2), yi);
dy = dy - round(dy./box(2))*box(2);
dy2 = (dy./sigma1).^2;
dz = bsxfun(@minus, data(:, 3), zi);
dz = dz - round(dz./box(3))*box(3);
dz2 = (dz./sigma1).^2;

gaussx = exp(-0.5 * dx2);
gaussy = exp(-0.5 * dy2);
gaussz = exp(-0.5 * dz2);

force = zeros(natom, 3);
t2 = zeros(nx, ny);
t3 = zeros(nx, ny, nz);

vec1 = zeros(nx, 1, 1);
vec2 = zeros(1, ny, 1);
vec3 = zeros(1, 1, nz);

for i = 1:natom_sub
  iatom = index_atom(i);
  t2 = bsxfun(@times, gaussx(i, :)', gaussy(i, :));
  vec3(1, 1, :) = gaussz(i, :);
  t3 = bsxfun(@times, t2, vec3);
  t3 = t3.*potential;

  vec1(:, 1, 1) = dx(i, :);
  vec2(1, :, 1) = dy(i, :);
  vec3(1, 1, :) = dz(i, :);

  tmp = bsxfun(@times, vec1, t3);
  force(iatom, 1) = sum(tmp(:));
  tmp = bsxfun(@times, vec2, t3);
  force(iatom, 2) = sum(tmp(:));
  tmp = bsxfun(@times, vec3, t3);
  force(iatom, 3) = sum(tmp(:));

  force(iatom, :) = force(iatom, :) * charge(iatom);
end

force = force * abs((xi(2)-xi(1)) * (yi(2)-yi(1)) * (zi(2)-zi(1))) / (2*pi*sqrt(2*pi)*(sigma1.^5));

