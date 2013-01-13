function rc = writenamdbin(filename, crd)
%% writenamdbin
% write namd restart (namdbin) file
%
%% Syntax
%# writenamdbin(filename, crd);
%
%% Description
% write namd restart (namdbin) file
% either coordinates or velocities
% are given
%
% * filename  - filename of namd restart (namdbin) file
% * crd       - coordinates [double 1 x natom3]
% * vel       - velocities [double 1 x natom3]
%
%% Example
%# writenamdbin('run.restart.coor',crd);
%# vel = vel / 20.45482706 % convert from Angstrom/ps to NAMD internal unit
%# writenamdbin('run.restart.vel',vel);
%
%% See also
% readnamdbin
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
nAtom = length(crd)/3;

%% write
disp('use little-endian')
fid = fopen(filename, 'w', 'l');

% disp('use big-endian')
% fid = fopen(filename,'w', 'b');

% nAtom
fwrite(fid, nAtom, 'int32');

% coordinates or velocities
fwrite(fid, crd, 'double');

fclose(fid);

