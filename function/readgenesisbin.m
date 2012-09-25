function [crd,vel,box] = readgenesisbin(filename, nAtomLimit)
%% readgenesisbin
% read genesis binary restart file
%
% function [crd,vel,box] = readgenesisbin(filename, nAtomLimit)
%
% input: filename ファイル名
%
% output: crd (1 x nAtom * 3) トラジェクトリ each row containing coordinates in the order [x1 y1 z1 x2 y2 z2 ...]
%       : vel (1 x nAtom * 3) トラジェクトリ each row containing coordinates in the order [x1 y1 z1 x2 y2 z2 ...]
%       : box (1 x 3) box-size
% 
% example:
% [crd,nAtom] = readnamdbin('eq.rst');
% 

fid = fopen(filename, 'r', 'l');

fseek(fid, 4, 0);
title1 = fread(fid, 80, '*char')';
fseek(fid, 4, 0);

fseek(fid, 4, 0);
title2 = fread(fid, 80, '*char')';
fseek(fid, 4, 0);

disp(title1);
disp(title2);

fseek(fid, 4, 0);
nAtom = fread(fid, 1, 'int32');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
rstfile_type = fread(fid, 1, 'int32');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
iseed = fread(fid, 1, 'int32');
num_deg_freedom = fread(fid, 1, 'int32');
fseek(fid, 4, 0);

box = zeros(1, 3);
fseek(fid, 4, 0);
box(1) = fread(fid, 1, 'double');
box(2) = fread(fid, 1, 'double');
box(3) = fread(fid, 1, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
thermostat_friction = fread(fid, 1, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
barostat_friction = fread(fid, 3, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
x = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
y = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
z = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
vx = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
vy = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

fseek(fid, 4, 0);
vz = fread(fid, nAtom, 'double');
fseek(fid, 4, 0);

crd = reshape([x'; y'; z'], 1, nAtom * 3);
vel = reshape([vx'; vy'; vz'], 1, nAtom * 3);

fclose(fid);


