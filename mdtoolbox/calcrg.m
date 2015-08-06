function rg = calcrg(trj, mass)
%% calcrg
% calc radius of gyration
% 
% example: 
% parm = readparm('ak.parm');
% trj  = readnetcdf('ak.nc');
% index_noh = ~selectname(parm.atom_name, 'H*');
% mass = parm.mass(index_noh);
% trj  = trj(:, to3(index_noh));
% rg   = calcrg(trj, mass);
% 

%% setup
nframe = size(trj, 1);
natom3 = size(trj, 2);
natom  = natom3/3;

if ~exist('mass', 'var') || isempty(mass)
  mass = ones(1, natom);
else
  if iscolumn(mass)
    mass = mass';
  end
end

%% calculation
rg = zeros(nframe, 1);
[trj, com] = decenter(trj, [], mass);
totalMass  = sum(mass);

x = sum(bsxfun(@times, trj(:, 1:3:end).^2, mass), 2);
y = sum(bsxfun(@times, trj(:, 2:3:end).^2, mass), 2);
z = sum(bsxfun(@times, trj(:, 3:3:end).^2, mass), 2);
rg = sqrt((x + y + z)./totalMass);

