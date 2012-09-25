function rc = writenamdbin(filename, crd)
%% writenamdbin
% write namd restart (namdbin) file
%
% function rc = writenamdbin(filename, crd)
%
% input: filename ファイル名
%        crd (1 x nAtom * 3) トラジェクトリ each row containing coordinates in the order [x1 y1 z1 x2 y2 z2 ...]
%
% output: rc
%
% example:
% writenamdbin('run.restart.coor',crd);
% vel = vel / 20.45482706 % convert from Angstrom/ps to NAMD internal unit
% writenamdbin('run.restart.vel',vel);
%

nAtom = length(crd)/3;

disp('use little-endian')
fid = fopen(filename, 'w', 'l');

% disp('use big-endian')
% fid = fopen(filename,'w', 'b');

% nAtom
fwrite(fid, nAtom, 'int32');

% coordinates or velocities
fwrite(fid, crd, 'double');

fclose(fid);


