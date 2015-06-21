function writegro(filename, gro, trj)
%% writegro
% write gromacs gro (Gromos87 format) file
%
%% Syntax
%# writegro(filename, gro);
%# writegro(filename, gro, crd);
%
%% Description
% This routine writes gromacs gro (Gromos87 format) file.
% Infomration of the gro file is stored in a structure variable.
% If coordinates are given, coordinates of output gro file 
% are repaced with the given coordinates
%
% * filename  - output filename of gromacs gro file [char]
% * gro       - structure data
%     residue_number: residue number [natom x 1 integer]
%       residue_name: residue name [natom x 5 char]
%          atom_name: atom name [natom x 5 char]
%        atom_number: atom number [natom x 1 integer]
%           position: atom XYZ coordinates [natom x 3 double]
%           velocity: atom XYZ velocities [natom x 3 double]
% * crd,trj   - coordinates [1 x natom3 double]
%               if the gro file contains multiple structures
%               coordinates [nstep x natom3 double]
% * vel       - velocities [1 x natom3 double]
%               if the gro file contains multiple structures
%               velocities [nstep x natom3 double]
% * box       - box size [1 x 3 double]
%
%% Example
%# gro = readpsf('jac.gro');
%# gro.xyz(:, 1) = gro.xyz(:, 1) + 2.5; %translate in x-axis by 2.5 Angstrom
%# writegro('jac_x.gro', gro);
%
%% See also
% readgro
% 
%% References
% http://manual.gromacs.org/current/online/gro.html
%
%% TODO
% velocities
% time step in title
% multiple steps
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% preparation
natom = size(gro.atom_number, 1);

if ~exist('trj', 'var') || isempty(trj)
  trj = gro.position';
  trj = trj(:)';
end
nstep = size(trj, 1);
  
%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write file
fprintf(fid, 'REMARKS FILENAME=%s CREATED BY MATLAB\n', filename);
fprintf(fid, '%5d\n', natom);

istep = 1;
for iatom = 1:natom
  fprintf(fid, '%5d', gro.residue_number(iatom));
  fprintf(fid, '%-5s', gro.residue_name(iatom, :));
  fprintf(fid, '%5s', gro.atom_name(iatom, :));
  fprintf(fid, '%5d', gro.atom_number(iatom, :));
  fprintf(fid, '%8.3f', trj(istep, 3*(iatom-1)+1));
  fprintf(fid, '%8.3f', trj(istep, 3*(iatom-1)+2));
  fprintf(fid, '%8.3f', trj(istep, 3*(iatom-1)+3));
  fprintf(fid, '\n');
end
  
fprintf(fid, '%10.5f%10.5f%10.5f\n', gro.box(1), gro.box(2), gro.box(3));

