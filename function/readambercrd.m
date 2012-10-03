function [crd, box, vel, title] = readambercrd(filename, index)
%% readambercrd
% read amber coordinate/restart file
%
%% Syntax
%# crd = readambercrd(filename);
%# crd = readambercrd(filename, index_atom);
%# [crd, box] = readambercrd(filename, index_atom);
%# [crd, box, vel] = readambercrd(filename, index_atom);
%# [crd, ~, vel] = readambercrd(filename, index_atom);
%# [crd, box, vel, title] = readambercrd(filename, index_atom);
%
%% Description
% The XYZ coordinates or velocities of atoms are read into 
% 'crd' and 'vel' variables which have '3*natom' columns
% in order [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input amber coordinate/restart filename
% * index_atom - atom index or logical index specifying atoms to be read
% * crd        - coordinates [1 x natom3 double]
% * box        - size of the periodic box [nstep x 3 double]
% * vel        - velocities [1 x natom3 double]
% * title      - title characters [chars]
%
%% Example
%# crd = readambercrd('ak.crd');
%
%% See alo
% writeambercrd
%
%% References
% http://ambermd.org/formats.html#restart
%

%% initialization
crd = [];
box = [];
vel = [];
title = [];

if nargin > 1
  if islogical(index)
    index = find(index);
  end
end

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse
title = fgetl(fid);
for i = (numel(title)+1):80
  title = [title ' '];
end
title = title(1:80);
str = fgetl(fid);
%natom = sscanf(str, '%5d', 1);
natom = sscanf(str, '%d', 1);
natom3 = natom*3;

if nargin == 1
  index = 1:natom;
end
index3 = to3(index);

%crd = fscanf(fid, '%12f', natom3);
crd = textscan(fid, '%12.7f', natom3);
crd = cell2mat(crd)';
crd = crd(1, index3);

vel = textscan(fid, '%12.7f', natom3);
vel = cell2mat(vel)';

if numel(vel) == 0
  vel = [];
elseif (numel(vel) == 3) | (numel(vel) == 6)
  box = vel(1:3);
  vel = [];
elseif numel(vel) == natom3
  vel = vel(1, index3);
  box = textscan(fid, '%12.7f', 6);
  box = cell2mat(box)';
  if numel(box) >= 3
    box = box(1:3);
  else
    box = [];
  end
end

