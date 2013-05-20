function kappa2 = calcorientationfactor(trj, index_donor, index_acceptor);
%% calcorientationfactor
% calculate the square of the orientation factor between two transition dipole moments
%
%% Syntax
%# kappa2 = calcorientationfactor(x, index_donor, index_acceptor);
%
%% Description
% square of the orientation factor (kappa2) has values between 0
% and 4. 
%
% * trj            - coordinates or trajectory
% * index_donor    - index of two donor-dipole atoms
% * index_acceptor - index of two acceptor-dipole atoms
% * kappa2         - square of orientation factor
%
%% Example
% trj = readdcd('run.dcd');
% kappa2 = calcorientationfactor(trj, [7 8], [20 21]);
%
%% See alo
% 

%% setup
nstep = size(trj, 1);
kappa = zeros(nstep, 1);

assert(numel(index_donor) == 2, 'index_donor must has 2 elements.')
assert(numel(index_acceptor) == 2, 'index_acceptor must has 2 elements.')

%% calculate orientation factor
% vector between point-like oscillating dipoles
[~, crd_donor] = decenter(trj(:, to3(index_donor)));
[~, crd_acceptor] = decenter(trj(:, to3(index_acceptor)));
r_da = crd_acceptor - crd_donor;
r_da = normalize_vector(r_da);

% dipole vectors
dipole_donor = trj(:, to3(index_donor(1))) - trj(:, to3(index_donor(2)));
dipole_donor = normalize_vector(dipole_donor);

dipole_acceptor = trj(:, to3(index_acceptor(1))) - trj(:, to3(index_acceptor(2)));
dipole_acceptor = normalize_vector(dipole_acceptor);

% angles
theta_donor  = dot(dipole_donor, r_da, 2);
theta_acceptor  = dot(dipole_acceptor, r_da, 2);
theta_da = dot(dipole_donor, dipole_acceptor, 2);

% orientation factor
kappa = theta_da - 3.0 .* theta_donor .* theta_acceptor;
kappa2 = kappa.^2;


%% subroutines
function y = normalize_vector(x)
x_amplitude = sqrt(dot(x, x, 2));
y = bsxfun(@rdivide, x, x_amplitude);


