function [crd, box, vel, title] = readmarblecrd(filename, index)
%% readmarblecrd
% read marble coordinate/restart file
%
%% Syntax
%# crd = readmarblecrd(filename);
%# crd = readmarblecrd(filename, index_atom);
%# [crd, box] = readmarblecrd(filename, index_atom);
%# [crd, box, vel] = readmarblecrd(filename, index_atom);
%# [crd, box, vel, title] = readmarblecrd(filename, index_atom);
%
%% Description
% The XYZ coordinates or velocities of atoms are read into 
% 'crd' and 'vel' variables which have '3*natom' columns
% in order [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input MARBLE coordinate/restart filename
% * index_atom - atom index or logical index specifying atoms to be read
% * crd        - coordinates [1 x natom3 double]
% * vel        - velocities [1 x natom3 double]
% * box        - size of the periodic box [nstep x 3 double]
% * title      - title characters [chars]
%
%% Example
%# crd = readmarblecrd('eq.crd');
%
%% See alo
% writemarbletrj
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
str = fgetl(fid);
%natom = sscanf(str, '%5d', 1);
natom = sscanf(str, '%d', 1);
natom3 = natom*3;

if nargin == 1
  index = 1:natom;
end
index3 = to3(index);

crd = textscan(fid, '%f', natom3);
crd = cell2mat(crd)';
crd = crd(1, index3);

vel = textscan(fid, '%f', natom3);
vel = cell2mat(vel)';

dummy = textscan(fid, '%d', 1);

box = textscan(fid, '%f', 9);
box = cell2mat(box)';
if numel(box) >= 9
  box = box([1 5 9]);
end

