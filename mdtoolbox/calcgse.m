function [energy, potential, force, xi, yi, zi] = calcgse(trj, charge, box, xi, yi, zi, sigma, weight);
%% calcgse
% compute smooth (or reciprocal) electrostatic potential, energy, and force using k-space Gaussian split Ewald
%
%% Syntax
%# [energy, potential, force, xi, yi, zi] = calcgse(trj, charge, box, xi, yi, zi, sigma, weight);
%
%% Description
%
%% Example
%# psf = readpsf('run.psf');
%# [trj, box] = readdcd('run.dcd');
%# xi = -32:32;
%# yi = -32:32;
%# zi = -32:32;
%# [energy, potential, force] = calcgse(trj, psf.charge, box, xi, yi, zi);
%# 
%
%% See alo
% 
% 

%% setup
coefficient = 332.0716;     % CHARMM
%coefficient = 332.05221729; % AMBER
%coefficient = 332.0636930;  % GROMACS

nstep  = size(trj, 1);
natom  = numel(charge);
energy = [];
force  = [];

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
  weight = ones(nstep, 1);
  weight = weight./sum(weight);
end


%% compute k-space Gaussian split Ewald
% pre-allocation
data = zeros(natom, 3);
potential = zeros(numel(yi), numel(xi), numel(zi));
energy = zeros(nstep, 1);
force = zeros(1, 3*natom);

for istep = 1:nstep
  % determine grids
  xi_gse = linspace(-box(istep, 1)/2, box(istep, 1)/2, floor(box(istep, 1))+1);
  yi_gse = linspace(-box(istep, 2)/2, box(istep, 2)/2, floor(box(istep, 2))+1);
  zi_gse = linspace(-box(istep, 3)/2, box(istep, 3)/2, floor(box(istep, 3))+1);
  xi_gse(end) = [];
  yi_gse(end) = [];
  zi_gse(end) = [];
  nx = numel(xi_gse);
  ny = numel(yi_gse);
  nz = numel(zi_gse);

  % 1: charge spreading by convolving Gaussian function in real-space
  data = reshape(trj(istep, :), 3, natom)';
  d1 = ksdensity3d(data, xi_gse, yi_gse, zi_gse, bandwidth1, box(istep, :), charge);
  d2 = permute(d1, [2,1,3]);

  % FFT
  f = fft3d(d2);

  % 2: on-mesh charge spreading by convolution with Gaussian function exp(-sig^2*k^2/2)
  [X, Y, Z] = meshgrid(2*pi*ifftshift(-ceil((nx-1)/2):floor((nx-1)/2))./box(istep, 1), ...
                     2*pi*ifftshift(-ceil((ny-1)/2):floor((ny-1)/2))./box(istep, 2), ...
                     2*pi*ifftshift(-ceil((nz-1)/2):floor((nz-1)/2))./box(istep, 3));
  gamma_mod = exp(-(sigma2.^2)*(X.^2 + Y.^2 + Z.^2)./2);

  % 3: solve poisson equation by convolution with Green function 4*pi/k^2
  gamma_mod = gamma_mod .* (4*pi) ./ (X.^2 + Y.^2 + Z.^2);
  gamma_mod(1, 1, 1) = 0;
  f = f.*gamma_mod;

  % 4: inverse FFT and get the on-mesh potential
  potential_gse = ifft3d(f);

  % 4': compute reciprocal energy
  d2 = d2.*potential_gse;
  energy(istep) = 0.5 * abs((xi_gse(2)-xi_gse(1)) * (yi_gse(2)-yi_gse(1)) * (zi_gse(2)-zi_gse(1))) * sum(d2(:));
  energy(istep) = energy(istep) * coefficient;

  % 4'': mesh interpolation if needed
  %potential = potential + weight*interp3(xi_gse, yi_gse, zi_gse, potential_gse, xi_query, yi_query, zi_query, 'linear');
  potential = potential + weight(istep)*interp3(xi_gse, yi_gse, zi_gse, potential_gse, xi_query, yi_query, zi_query, 'cubic');
end

potential = permute(potential, [2,1,3]);

%%%%%% 3D FFT routines
function A_k = fft3d(A_r);
A_k = fft(fft(fft(A_r,[],1),[],2),[],3);

function A_r = ifft3d(A_k);
A_r = ifft(ifft(ifft(A_k,[],1),[],2),[],3,'symmetric');

