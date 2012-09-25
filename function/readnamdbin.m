function [crd, nAtom] = readnamdbin(filename, nAtomLimit)
% read namd restart (namdbin) file
%
% function [crd,nAtom] = readnamdbin(filename)
%
% input: filename ファイル名
%
% output: crd (1 x nAtom * 3) トラジェクトリ each row containing coordinates in the order [x1 y1 z1 x2 y2 z2 ...]
%       : nAtom 読み込まれた原子数
% 
% example:
% [crd,nAtom] = readnamdbin('run.restart.coor');
% [vel,nAtom] = readnamdbin('run.restart.vel');
% vel = vel * 20.45482706 % convert to Angstrom/ps
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


