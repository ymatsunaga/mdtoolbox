function writenmd(filename, crd, mode, resid, chainid)
%% writenmd
% write nmd format file for Normal Mode Wizard in VMD
%
%% Syntax
%# writenmd(filename, crd, mode);
%# writenmd(filename, crd, mode, resid);
%# writenmd(filename, crd, mode, resid, chainid);
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
% * resids    - residue numbers [optional] [nstep3 x m char]
% * chainids  - chain identifiers [optional] [nstep3 x 1 integer]
%
%% Example
%# [pdb, crd] = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb = substruct(pdb, index_ca);
%# crd = crd(to3(index_ca));
%# mode = anm(crd, 8.0);
%# writenmd('anm.nmd', crd, mode)
%# In VMD
%# nmwiz load anm.nmd
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
is_resid   = false;
is_chainid = false;

natom3 = numel(crd);
natom  = natom3/3;
nmode = size(mode, 2);

if exist('resid', 'var')
  is_resid = true;
end

if exist('chainid', 'var')
  is_chainid = true;
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

if is_chainid
  fprintf(fid, 'chids');
  for i = 1:natom
    fprintf(fid, ' %s', chainid(i, :));
  end
  fprintf(fid, '\n');
end

if is_resid
  fprintf(fid, 'resnums');
  for i = 1:natom
    fprintf(fid, ' %d', resid(i));
  end
  fprintf(fid, '\n');
end

for imode = 1:nmode
  fprintf(fid, 'mode');
  fprintf(fid, ' %d', imode);
  %fprintf(fid, ' %f', norm(mode(:, imode)));
  fprintf(fid, ' %f', 1.0);
  for i = 1:natom3
    fprintf(fid, ' %f', mode(i, imode));
  end
  fprintf(fid, '\n');
end

