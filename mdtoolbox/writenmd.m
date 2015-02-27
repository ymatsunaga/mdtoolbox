function writenmd(filename, crd, mode, nmode)
%% writenmd
% write nmd format file for Normal Mode Wizard in VMD
%
%% Syntax
%# writenmd(filename, crd, mode)
%# writenmd(filename, crd, mode, nmode)
%
%% Description
% Output NMD file for visualization of Normal Mode
% by using Normal Mode Wizard in VMD. 
% For details about the Normal Mode Wizard, 
% see http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/
%
% * filename  - output nmd filename
% * crd       - coordinates [1 x natom3 double]
% * mode      - mode information [nstep3 x n double]
% * nmode     - number of modes to be written [integer]
%               if omitted, all modes in 'mode' is written.
%
%% Example
%# [pdb, crd] = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb = substruct(pdb, index_ca);
%# crd = crd(to3(index_ca));
%# crd = decenter(crd);
%# emode = anm(crd, 12.0);
%# writenmd(filename, crd, mode, 10)
%
%% See also
% anm, calcpca
%
%% References
% http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/#nmd
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% setup
natom3 = numel(crd);
natom  = natom3/3;

if ~exist('nmode', 'var') || isempty(nmode)
  nmode = size(mode, 2);
end

%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write data
fprintf(fid, 'coordinates');
for i = 1:natom3
  fprintf(fid, ' %f', crd(i));
end
fprintf(fid, '\n');

for imode = 1:nmode
  fprintf(fid, 'mode');
  fprintf(fid, ' %d', imode);
  fprintf(fid, ' 1.0');
  for i = 1:natom3
    fprintf(fid, ' %f', mode(i, imode));
  end
  fprintf(fid, '\n');
end

