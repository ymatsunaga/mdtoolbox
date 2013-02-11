function mass = calcmass(name);
%% calcmass
% calculate to the masses of input atom names
%
%% Syntax
%# mass = calcmass(name);
%
%% Description
% Calculate to the masses of input atom names
%
% * name    - atom names [natom x n chars]
% * mass    - masses [natom x 1 double]
%
%% Example
%# parm = readparm('ak.parm');
%# mass = calcmass(parm.atom_name)
%
%% See alo
% 

%% setup
natom       = size(name,1);
mass        = zeros(natom, 1);
index_check = false(natom, 1);

%% convert to mass
% basic atoms
index = selectname(name, 'H*');
index_check = index_check | index;
mass(index) = 1.00800;

index = selectname(name, 'C*');
index_check = index_check | index;
mass(index) = 12.01100;

index = selectname(name, 'N*');
index_check = index_check | index;
mass(index) = 14.00700;

index = selectname(name, 'O*');
index_check = index_check | index;
mass(index) = 15.99900;

index = selectname(name, 'S*');
index_check = index_check | index;
mass(index) = 32.06000;

% ions, etc.
index = selectname(name, 'LIT');
index_check = index_check | index;
mass(index) = 6.94100;

index = selectname(name, 'SOD');
index_check = index_check | index;
mass(index) = 22.98977;

index = selectname(name, 'MG');
index_check = index_check | index;
mass(index) = 24.30500;

index = selectname(name, 'POT');
index_check = index_check | index;
mass(index) = 39.09830;

index = selectname(name, 'CAL');
index_check = index_check | index;
mass(index) = 40.08000;

index = selectname(name, 'RUB');
index_check = index_check | index;
mass(index) = 85.46780;

index = selectname(name, 'CES');
index_check = index_check | index;
mass(index) = 132.90545;

index = selectname(name, 'BAR');
index_check = index_check | index;
mass(index) = 137.32700;

index = selectname(name, 'ZN');
index_check = index_check | index;
mass(index) = 65.37000;

index = selectname(name, 'CAD');
index_check = index_check | index;
mass(index) = 112.41100;

index = selectname(name, 'CLA');
index_check = index_check | index;
mass(index) = 35.45000;

index = selectname(name, 'OT');
index_check = index_check | index;
mass(index) = 15.99940;

index = selectname(name, 'OX');
index_check = index_check | index;
mass(index) = 15.99940;

index = selectname(name, 'PL');
index_check = index_check | index;
mass(index) = 30.974000;

index = selectname(name, 'SL');
index_check = index_check | index;
mass(index) = 32.060000;

if any(~index_check)
  fprintf('could not determine the masses of following atoms:\n');
  disp(name(index_check == false, :));
end

