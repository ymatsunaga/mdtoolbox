function [crd, vel, box, header] = readgenesisbin(filename)
%% readgenesisbin
% read genesis binary restart file
%
%% Syntax
%# crd = readgenesisbin(filename);
%# [crd, vel] = readgenesisbin(filename);
%# [crd, vel, box] = readgenesisbin(filename);
%
%% Description
% read genesis binary restart file
%
% * filename  - filename of genesis biary restart
% * crd       - coordinates [1 x natom3]
% * vel       - velocities [1 x natom3]
% * box       - box size [1 x 3]
%
%% Example
%# crd = readgenesisbin('run.rst');
%
%% See also
% writegenesisbin
%

%% check endian
% try little endian
assert(ischar(filename), 'Please specify valid filename for the first argument')
fid = fopen(filename, 'r', 'l');
ival = fread(fid, 1, 'int32');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% try big endian
if ival ~= 80
  fclose(fid);
  fid = fopen(filename, 'r', 'l');
  ival = fread(fid, 1, 'int32');
end

% try little-endian 64-bits
if ival ~= 80
  fclose(fid);
  fid = fopen(filename, 'r', 'a');
  ival = fread(fid, 1, 'int32');
end

% try big-endian 64-bits
if ival ~= 80
  fclose(fid);
  fid = fopen(filename, 'r', 's');
  ival = fread(fid, 1, 'int32');
end

%% read data
fseek(fid, 0, 'bof');

fseek(fid, 4, 0);
header.title1 = fread(fid, 80, '*char')';
fseek(fid, 4, 0);

fseek(fid, 4, 0);
header.title2 = fread(fid, 80, '*char')';
fseek(fid, 4, 0);

disp(header.title1);
disp(header.title2);

fseek(fid, 4, 0);
header.nAtom = fread(fid, 1, 'int32');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
header.rstfile_type = fread(fid, 1, 'int32');
fseek(fid, 4, 0);

%% type == minimization
if header.rstfile_type == 1
  box = zeros(1, 3);
  fseek(fid, 4, 0);
  box(1) = fread(fid, 1, 'double');
  box(2) = fread(fid, 1, 'double');
  box(3) = fread(fid, 1, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  header.energy = fread(fid, 1, 'double');
  header.delta  = fread(fid, 1, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  x = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  y = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  z = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  crd = reshape([x'; y'; z'], 1, header.nAtom * 3);

%% type == MD
elseif header.rstfile_type == 2
  fseek(fid, 4, 0);
  header.iseed = fread(fid, 1, 'int32');
  header.num_deg_freedom = fread(fid, 1, 'int32');
  fseek(fid, 4, 0);

  box = zeros(1, 3);
  fseek(fid, 4, 0);
  box(1) = fread(fid, 1, 'double');
  box(2) = fread(fid, 1, 'double');
  box(3) = fread(fid, 1, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  header.thermostat_momentum = fread(fid, 1, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  header.barostat_momentum = fread(fid, 3, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  x = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  y = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  z = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  vx = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  vy = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);

  fseek(fid, 4, 0);
  vz = fread(fid, header.nAtom, 'double');
  fseek(fid, 4, 0);
  
  fseek(fid, 4, 0);
  header.nrandom = fread(fid, 1, 'int32');
  fseek(fid, 4, 0);

  if header.nrandom > 0
    fseek(fid, 4, 0);
    header.random = fread(fid, header.nrandom, 'double');
    fseek(fid, 4, 0);
  end

  crd = reshape([x'; y'; z'], 1, header.nAtom * 3);
  vel = reshape([vx'; vy'; vz'], 1, header.nAtom * 3);

%% type == REMD
elseif header.rstfile_type == 3
  error('REMD type is not supported');
  
else
  error(sprintf('%d is undefined type', header.rstfile_type));

end



