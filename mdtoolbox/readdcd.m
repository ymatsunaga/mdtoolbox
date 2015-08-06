function [trj, box, header] = readdcd(filename, index)
%% readdcd
% read xplor or charmm (namd) format dcd file
%
%% Syntax
%# trj = readdcd(filename);
%# trj = readdcd(filename, index_atom);
%# [trj, box] = readdcd(filename, index_atom);
%# [trj, box, header] = readdcd(filename, index_atom);
%# [trj, ~, header] = readdcd(filename, index_atom);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nframe' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input dcd trajectory filename
% * index_atom - index or logical index for specifying atoms to be read
% * trj        - trajectory [nframe x natom3 double]
% * box        - box size [nframe x 3 double]
% * header     - structure variable, which has header information 
%                [structure]
%
%% Example
%# trj = readdcd('ak.dcd');
%
%% See also
% writedcd
%
%% References for dcd format
% MolFile Plugin http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
% CafeMol Manual http://www.cafemol.org/doc.php
% EGO_VIII Manual http://www.lrz.de/~heller/ego/manual/node93.html
%
%% TODO
% try-catch
%

%% initialization
trj = [];
box = [];
header.is_charmm = false;
header.is_charmm_extrablock = false;
header.is_charmm_4dims = false;

if exist('index', 'var') && ~isempty(index)
  if islogical(index)
    index = find(index);
  end
end

%% open file
assert(ischar(filename), 'Please specify valid filename for the first argument');
fid = fopen(filename, 'r', 'l'); % try little-endian
assert(fid > 0, 'Could not open file.');

% check the size of file
fseek(fid, 0, 'eof');
filesize = ftell(fid); % get file size in byte
fseek(fid, 0, 'bof');

%% read block 1 (header)
% read blocksize1 and check endian (seek=0 byte)
header.blocksize1 = fread(fid, 1, 'int32');
if header.blocksize1 ~= 84
  fclose(fid);
  fid = fopen(filename, 'r', 'b'); % try big-endian
  header.blocksize1 = fread(fid, 1, 'int32');
end

if header.blocksize1 ~= 84
  fclose(fid);
  fid = fopen(filename, 'r', 'a'); % try little-endian 64-bits
  header.blocksize1 = fread(fid, 1, 'int32');
end

if header.blocksize1 ~= 84
  fclose(fid);
  fid = fopen(filename, 'r', 's'); % try big-endian 64-bits
  header.blocksize1 = fread(fid, 1, 'int32');
end

if header.blocksize1 ~= 84
  error('file may not be a dcd file.');
end

cleaner = onCleanup(@() fclose(fid));

% header (4 chars) either "CORD" or "VELD" (seek=4 byte)
header.hdr = fread(fid, 4, 'char');
header.hdr = char(header.hdr)';

% the total # of frames (snapshots) (seek=8 byte)
header.nset = fread(fid, 1, 'int32');

% starting time-step (seek=12 byte)
header.istrt = fread(fid, 1, 'int32');

% frequency to save trajectory (seek=16 byte)
header.nsavc = fread(fid, 1, 'int32');

% the total # of simulation step (seek=20 byte)
header.nstep = fread(fid, 1, 'int32');

% null4 (int*4) (seek=24 byte)
header.null4 = fread(fid, 4, 'int32');

% # of degrees of freedom (seek=40 byte)
header.nfreat = fread(fid, 1, 'int32');

% step-size of simulation (seek=44 byte)
header.delta = fread(fid, 1, 'float32');

% null9 (int*9) (seek=48 byte)
header.null9 = fread(fid, 9, 'int32');

% version (seek=84 byte)
header.version = fread(fid, 1, 'int32');

% charmm extension format
if header.version > 0
  header.is_charmm = true;
end

% delta is double precision in xplor format
if ~header.is_charmm
  cof = ftell(fid);

  % reread delta in double precision
  fseek(fid, 44, -1);
  header.delta = fread(fid, 1, 'float64');
  
  fseek(fid, cof, -1);
end  

% check charmm extensions
if header.is_charmm
  cof = ftell(fid);

  % charmm extrablock extension
  fseek(fid, 48, -1);
  n = fread(fid, 1, 'int32');
  if n == 1
    header.is_charmm_extrablock = true;
  end

  % charmm 4dims extension
  fseek(fid, 52, -1);
  n = fread(fid, 1, 'int32');
  if n == 1
    header.is_charmm_extrablock = true;
  end

  fseek(fid, cof, -1);
end  

% blocksize1
blocksize1 = fread(fid, 1, 'int32');
assert(header.blocksize1 == blocksize1, 'blocksize1 values are not consistent');

%% read block 2 (title)
% blocksize2
header.blocksize2 = fread(fid, 1, 'int32');

% # of title lines
header.ntitle = fread(fid, 1, 'int32');

% title
header.title = fread(fid, 80*header.ntitle, 'char');
header.title = char(reshape(header.title, 80, [])');

% blocksize2
blocksize2 = fread(fid, 1, 'int32');
assert(header.blocksize2 == blocksize2, 'blocksize2 values are not consistent');

%% read block 3 (natom)
% blocksize3
header.blocksize3 = fread(fid, 1, 'int32');

% # of atoms
header.natom = fread(fid, 1, 'int32');

% blocksize3
blocksize3 = fread(fid, 1, 'int32');
assert(header.blocksize3 == blocksize3, 'blocksize3 values are not consistent');

%% read coordinates
headersize = ftell(fid);

if header.is_charmm_extrablock
  extrablocksize = 4*2 + 8*6;
else
  extrablocksize = 0;
end

coordblocksize = (4*2 + 4*header.natom)*3;

nframe = floor(filesize - headersize) / (extrablocksize + coordblocksize);

if ~exist('index', 'var') || isempty(index)
  index = 1:header.natom;
end

trj = zeros(nframe, numel(index)*3);
box = zeros(nframe, 3);

% read next frames
for iframe = 1:nframe
  % read charmm extrablock (unitcell info)
  if header.is_charmm_extrablock
    blocksize = fread(fid, 1, 'int32');
    dummy = fread(fid, 6, 'float64');
    blocksize = fread(fid, 1, 'int32');
  end

  % read x coordinates
  blocksize = fread(fid, 1, 'int32');
  x = fread(fid, blocksize/4, 'float32');
  blocksize = fread(fid, 1, 'int32');

  % read y coordinates 
  blocksize = fread(fid, 1, 'int32');
  y = fread(fid, blocksize/4, 'float32');
  blocksize = fread(fid, 1, 'int32');

  % read z coordinates 
  blocksize = fread(fid, 1, 'int32');
  z = fread(fid, blocksize/4, 'float32');
  blocksize = fread(fid, 1, 'int32');
  
  % ignore charmm 4dims extension
  if header.is_charmm_4dims
    blocksize = fread(fid, 1, 'int32');
    fseek(fid, blocksize, 0);
    blocksize = fread(fid, 1, 'int32');
  end

  if header.is_charmm_extrablock
    box(iframe, :) = dummy([1 3 6])';
  end
  trj(iframe, 1:3:end) = x(index)';
  trj(iframe, 2:3:end) = y(index)';
  trj(iframe, 3:3:end) = z(index)';
end

