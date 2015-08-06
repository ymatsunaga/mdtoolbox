function [gro, crd, vel] = readgro(filename, index)
%% readgro
% read gromacs gro (Gromos87 format) file
%
%% Syntax
%# [gro, crd] = readgro(filename);
%# [gro, crd, vel] = readgro(filename);
%
%% Description
% This routine reads gromacs gro (Gromos87 format) file.
% Infomration of the gro file is stored in a structure variable.
% Coordinates, and velocities can be separately obtained
% if the user specifies additional output arguments.
%
% * filename  - input filename of gromacs gro file [char]
% * gro       - structure data
%     residue_number: residue number [natom x 1 integer]
%       residue_name: residue name [natom x 5 char]
%          atom_name: atom name [natom x 5 char]
%        atom_number: atom number [natom x 1 integer]
%           position: atom XYZ coordinates [natom x 3 double]
%           velocity: atom XYZ velocities [natom x 3 double]
% * crd       - coordinates [1 x natom3 double]
%               if the gro file contains multiple structures
%               coordinates [nframe x natom3 double]
% * vel       - velocities [1 x natom3 double]
%               if the gro file contains multiple structures
%               velocities [nframe x natom3 double]
% * box       - box size [1 x 3 double]
%
%% Example
%# [gro, crd] = readgro('out.gro');
%
%% See alo
% writegro
%
%% References
% http://manual.gromacs.org/current/online/gro.html
%
%% TODO
% velocities
% time step in title
% multiple frames
%

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse
title = fgetl(fid);
str = fgetl(fid);
natom = sscanf(str, '%d', 1);

%% allocate variables
gro.residue_number = zeros(natom, 1);
gro.residue_name   = repmat('12345', natom, 1);
gro.atom_name      = repmat('12345', natom, 1);
gro.atom_number    = zeros(natom, 1);
gro.position       = zeros(natom, 3);
gro.velocity       = zeros(natom, 3);
gro.box            = zeros(1, 3);

%% read file
for iatom = 1:natom
  line = fgetl(fid);

  gro.residue_number(iatom)  = str2num(line(1:5));
  gro.residue_name(iatom, :) = line(6:10);
  gro.atom_name(iatom, :)    = line(11:15);
  gro.atom_number(iatom)     = str2num(line(16:20));
  gro.position(iatom, 1)     = str2num(line(21:28));
  gro.position(iatom, 2)     = str2num(line(29:36));
  gro.position(iatom, 3)     = str2num(line(37:44));
  
end

line = fgetl(fid);
gro.box(1, :) = sscanf(line, '%f', 3);

if nargout >= 2
  crd = zeros(1, 3*natom);
  tmp = gro.position';
  crd(1, :) = tmp(:)';
end

