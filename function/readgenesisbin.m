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

fid = fopen(filename, 'r', 'l');

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
header.thermostat_friction = fread(fid, 1, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
header.barostat_friction = fread(fid, 3, 'double');
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

crd = reshape([x'; y'; z'], 1, header.nAtom * 3);
vel = reshape([vx'; vy'; vz'], 1, header.nAtom * 3);

fclose(fid);


