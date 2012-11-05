function [crd, nAtom] = readnamdbin(filename, nAtomLimit)
%% readnamdbin
% read namd restart (namdbin) file
%
%% Syntax
%# crd = writenamdbin(filename);
%# [crd, natom] = writenamdbin(filename);
%# [crd, natom] = writenamdbin(filename, natomlimit);
%
%% Description
% read namd restart (namdbin) file
%
% * filename   - filename of namd restart (namdbin) file
% * crd        - coordinates or velocities [double 1 x natom3]
% * natomlimit - maximum number of atoms to be read [integer]
% * natom      - number of atoms [integer]
%
%
%% Example
%# crd = readnamdbin('run.restart.coor');
%# vel = readnamdbin('run.restart.vel');
%# vel = vel * 20.45482706 % convert to Angstrom/ps
%
%% See also
% writenamdbin
%

if nargin == 1
  nAtomLimit = 10^9;
end

fid = fopen(filename,'r', 'l');

% check endian and read nAtom
iCheckEndian = fread(fid, 1, 'int32');
if iCheckEndian < nAtomLimit
  disp('assume big-endian');
else
  disp('assume little-endian')
  fclose(fid);
  fid = fopen(filename, 'r', 'l');
  iCheckEndian = fread(fid, 1, 'int32');
  if iCheckEndian > nAtomLimit
    error('can not tell the endianism...')
  end
end

nAtom = iCheckEndian;

% coordinates or velocities
crd = fread(fid, nAtom*3, 'double');

fclose(fid);


