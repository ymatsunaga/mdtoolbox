function [potential, xi, yi, zi, energy, force] = calcgse(trj, charge, box, sigma, weight);
%% calcgse
% compute smoothed (or reciprocal) electrostatic interaction energy, potential, and force using k-space Gaussian split Ewald
%
%% Syntax
%# [potential, xi, yi, zi, energy, force] = calcgse(trj, charge, box, sigma, weight);
%
%% Description
%
%% Example
%# psf = readpsf('run.psf');
%# [trj, box] = readdcd('run.dcd');
%# [potential, xi, yi, zi, energy] = calcgse(trj, psf.charge, box);
%# 
%
%% See alo
% 
% 

%% setup
coefficient = 332.0716;
nstep  = size(trj, 1);
natom  = numel(charge);
energy = [];
force = [];

if ~exist('sigma', 'var') || isempty(sigma)
  alpha = 0.34;
  sigma = 1./(alpha*sqrt(2));
  fprintf('message: sigma of %f A is used\n', sigma);
end
kappa = 0.35;
%kappa = 0.45;
%kappa = 0.5; % 0 <= kappa <= 0.5
sigma1 = sigma*sqrt(kappa);
sigma2 = sigma*sqrt(1 - 2*kappa);
bandwidth1 = [sigma1 sigma1 sigma1];

istep = 1;

%% determine grids
xi = linspace(-box(istep, 1)/2, box(istep, 1)/2, floor(box(istep, 1))+1);
yi = linspace(-box(istep, 2)/2, box(istep, 2)/2, floor(box(istep, 2))+1);
zi = linspace(-box(istep, 3)/2, box(istep, 3)/2, floor(box(istep, 3))+1);
xi(end) = [];
yi(end) = [];
zi(end) = [];
nx = numel(xi);
ny = numel(yi);
nz = numel(zi);
fprintf('message: %d grids in x-axis\n', nx);
fprintf('message: %d grids in y-axis\n', ny);
fprintf('message: %d grids in z-axis\n', nz);

%% 1: charge spreading by convolving Gaussian function in real-space
data = reshape(trj(istep, :), 3, natom)';
d1 = ksdensity3d(data, xi, yi, zi, bandwidth1, box(istep, :), charge);

%% FFT
f = fft3d(d1);

%% 2: on-mesh charge spreading by convolution with Gaussian function exp(-sig^2*k^2/2)
[X, Y, Z] = meshgrid(2*pi*ifftshift(-ceil((ny-1)/2):floor((ny-1)/2))./box(istep, 2), ...
                     2*pi*ifftshift(-ceil((nx-1)/2):floor((nx-1)/2))./box(istep, 1), ...
                     2*pi*ifftshift(-ceil((nz-1)/2):floor((nz-1)/2))./box(istep, 3));
gamma_mod = exp(-(sigma2.^2)*(X.^2 + Y.^2 + Z.^2)./2);

%% 3: solve poisson equation by convolution with Green function 4*pi/k^2
gamma_mod = gamma_mod .* (4*pi) ./ (X.^2 + Y.^2 + Z.^2);
gamma_mod(1, 1, 1) = 0;
f = f.*gamma_mod;

%% 4: inverse FFT and get the on-mesh potential
potential = ifft3d(f);

%% 4': compute reciprocal energy
d1 = d1.*potential;
energy = 0.5 * abs(xi(2)-xi(1)) * abs(yi(2)-yi(1)) * abs(zi(2)-zi(1)) * sum(d1(:));
energy = energy * coefficient;


%%%%%% 3D FFT routines
function A_k = fft3d(A_r);
A_k = fft(fft(fft(A_r,[],1),[],2),[],3);

function A_r = ifft3d(A_k);
A_r = ifft(ifft(ifft(A_k,[],1),[],2),[],3,'symmetric');
%A_r = ifft(ifft(ifft(A_k,[],1),[],2),[],3);

