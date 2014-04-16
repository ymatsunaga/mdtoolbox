function writedcd(filename, trj, box, header)
%% writedcd
% write xplor or charmm (namd) format dcd file
%
%% Syntax
%# writedcd(filename, trj);
%# writedcd(filename, trj, box);
%# writedcd(filename, trj, box, header);
%# writedcd(filename, trj, [], header);
%
%% Description
% This code puts trajectories into a dcd file. 
% If header information is not given, 
% default values are assumed. 
%
% * filename  - output dcd trajectory filename
% * trj       - trajectory [nstep x natom3 double]
% * box       - box size [nstep x 3 double]
% * header    - structure variable, which has header information 
%               [structure]
%
%% Example
%# trj = readdcd('ak.dcd');
%# trj = trj(:, 1:3:end) + 1.5;
%# writedcd('ak_translated.dcd', trj);
%
%% See also
% readdcd
%
%% References
% MolFile Plugin http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
% CafeMol Manual http://www.cafemol.org/doc.php
% EGO_VIII Manual http://www.lrz.de/~heller/ego/manual/node93.html
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
[nstep, natom3] = size(trj);
natom = natom3 / 3;

if ~exist('header', 'var') || isempty(header)
  % default header in xplor format
  header.is_charmm = false;
  header.is_charmm_extrablock = false;
  header.is_charmm_4dims = false;
  header.blocksize1 = 84;
  header.hdr = 'CORD';
  header.nset = size(trj, 1);
  header.istrt = 0;
  header.nsavc = 1;
  header.nstep = 0;
  header.null4 = zeros(4, 1);
  header.nfreat = 0;
  header.delta = 1.0;
  header.null9 = zeros(9, 1);
  header.version = 0;
  header.blocksize2 = 164;
  header.ntitle = 2;
  title1 = sprintf('REMARKS FILENAME=%s CREATED BY MATLAB', filename);
  for i = (numel(title1)+1):80
    title1 = [title1 ' '];
  end
  title1 = title1(1:80);
  title2 = sprintf('REMARKS DATE: %s CREATED BY USER: %s', datestr(now, 'mm/dd/yy'), getenv('USER'));
  for i = (numel(title2)+1):80
    title2 = [title2 ' '];
  end
  title2 = title2(1:80);
  header.title = [title1; title2];
  header.blocksize3 = 4;
  header.natom = size(trj, 2) / 3;
end

if header.nset ~= size(trj, 1)
  header.nset = size(trj, 1);
end

if exist('box', 'var') && ~isempty(box)
  % charmm format
  header.is_charmm = true;
  header.is_charmm_extrablock = true;
  if header.version == 0
    header.version = 1; % is_charmm -> true
  end
  header.null9(1) = 1; % is_charmm_extrablock -> true
else
  % xplor format
  header.is_charmm = false;
  header.is_charmm_extrablock = false;
  header.is_charmm_4dims = false;
  header.version = 0; % is_charmm -> false
  header.null9(1) = 0; % is_charmm_extrablock -> false
  header.null9(2) = 0; % is_charmm_4dims -> false
end

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write block 1 (header)
fwrite(fid, header.blocksize1, 'int32');             
fwrite(fid, header.hdr, 'uchar');
fwrite(fid, header.nset, 'int32');
fwrite(fid, header.istrt, 'int32');
fwrite(fid, header.nsavc, 'int32');
fwrite(fid, header.nstep, 'int32');
fwrite(fid, header.null4, 'int32');
fwrite(fid, header.nfreat, 'int32');

if header.is_charmm
  % charmm format
  fwrite(fid, header.delta, 'float32');
  fwrite(fid, header.null9, 'int32');
else
  % xplor format
  fwrite(fid, header.delta, 'float64');
  fwrite(fid, header.null9(2:end), 'int32');
end

fwrite(fid, header.version, 'int32');
fwrite(fid, header.blocksize1, 'int32');

%% write block 2 (title)
fwrite(fid, header.blocksize2, 'int32'); 
fwrite(fid, header.ntitle, 'int32');
fwrite(fid, header.title(1, :), 'uchar');
fwrite(fid, header.title(2, :), 'uchar');
fwrite(fid, header.blocksize2, 'int32');

%% write block 3 (natom)
fwrite(fid, header.blocksize3, 'int32');
fwrite(fid, header.natom, 'int32');
fwrite(fid, header.blocksize3, 'int32');

%% write coordinates
dummy = zeros(1, 6);
for istep = 1:nstep
  if header.is_charmm_extrablock
    fwrite(fid, 48, 'int32');
    dummy(1, [1 3 6]) = box(istep, :);
    fwrite(fid, dummy, 'float64');
    fwrite(fid, 48, 'int32');
  end
  
  fwrite(fid, natom*4, 'int32');
  fwrite(fid, trj(istep, 1:3:end), 'float32');
  fwrite(fid, natom*4, 'int32');

  fwrite(fid, natom*4, 'int32');
  fwrite(fid, trj(istep, 2:3:end), 'float32');
  fwrite(fid, natom*4, 'int32');

  fwrite(fid, natom*4, 'int32');
  fwrite(fid, trj(istep, 3:3:end), 'float32');
  fwrite(fid, natom*4, 'int32');
end

