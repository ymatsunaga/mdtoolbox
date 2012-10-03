function rc = writenamdbin(filename, crd)
%% writenamdbin
% write namd restart (namdbin) file
%
%% Syntax
%# readgenesisbin(filename, crd);
%# readgenesisbin(filename, vel);
%
%% Description
% write namd restart (namdbin) file
% either coordinates or velocities
% are given
%
% * filename  - filename of namd restart (namdbin) file
% * crd       - coordinates [1 x natom3]
% * vel       - velocities [1 x natom3]
%
%% Example
%# writenamdbin('run.restart.coor',crd);
%# vel = vel / 20.45482706 % convert from Angstrom/ps to NAMD internal unit
%# writenamdbin('run.restart.vel',vel);
%
%% See also
% readnamdbin
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

