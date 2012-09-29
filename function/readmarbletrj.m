function [trj, box, vel] = readmarbletrj(filename,index)
%% readmarbletrj
% read marble ascii-format trajectory file
%
%% Syntax
%# trj = readmarbletrj(filename);
%# trj = readmarbletrj(filename, index);
%# [trj, box] = readmarbletrj(filename, index);
%# [trj, box, vel] = readmarbletrj(filename, index);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nstep' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename  - input marble trajectory filename
% * index     - index or logical index specifying atoms to be read
% * trj       - trajectory [nstepx3natom double]
% * box       - box size [nstepx3 double]
% * vel       - velocities [nstepx3natom double]
%
%% Example
%# trj = readmarbletrj('eq.trj');
%
%% See alo
% readmarblecrd
%

%% initialization
trj = [];
box = [];
vel = [];
is_trj = false;
is_box = false;
is_vel = false;
iblock = 1;
is_compressed = false;

%% open file
filename = strtrim(filename);
if (numel(filename) >= 3) & strncmpi(filename((end-2):end), '.gz', numel('.gz'))
  dirname = tempname();
  dirname = [dirname '/'];
  mkdir(dirname);
  disp(sprintf('uncompressing %s into %s', filename, dirname))
  filename = gunzip(filename, dirname);
  filename = filename{1};
  disp('done')
  cleaner_rmdir = onCleanup(@() rmdir(dirname, 's'));
end

fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% read header
line = strtrim(fgetl(fid));
line = regexp(line, ' ', 'split');
natom = str2num(line{6}(1:end-1));
natom3 = natom*3;
trj_type = line{8};

if regexp(trj_type, 'X')
  is_trj = true;
end

if regexp(trj_type, 'B')
  is_box = true;
end

if regexp(trj_type, 'V')
  is_vel = true;
end

if nargin < 2
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
end
index3 = to3(index);

% buffer size is about 1 GByte
nblock = ceil(10^9 / (8*numel(index3)));
if nblock < 3
  nblock = 3;
end

% allocate buffer
if is_trj
  trj_buffer = zeros(nblock, numel(index3));
end

if is_box
  box_buffer = zeros(nblock, 3);
end

if is_vel
  vel_buffer = zeros(nblock, numel(index3));
end

%% parse
while ~feof(fid)
  if is_trj
    %crd = fscanf(fid, '%f %f %f\n', [3 natom]);
    %trj_buffer(iblock, :) = crd';
    crd = textscan(fid, '%f', natom3);
    xx = cell2mat(crd);
    if numel(xx) < natom3; break; end
    trj_buffer(iblock, :) = xx(index3)';
  end

  if is_box
    %crd2 = fscanf(fid, '%f %f %f\n', [3 3]);
    %box_buffer(iblock, :) = diag(crd2)';
    crd2 = textscan(fid, '%f', 9);
    bb = cell2mat(crd2);
    if numel(bb) < 9; break; end
    bb = bb([1 5 9]);
    box_buffer(iblock, :) = bb';
  end

  if is_vel
    crd = textscan(fid, '%f', natom3);
    xx = cell2mat(crd);
    if numel(xx) < natom3; break; end
    vel_buffer(iblock, :) = xx(index3)';
  end

  iblock = iblock + 1;
  
  if iblock > nblock
    if is_trj
      trj = [trj; trj_buffer];
    end
    if is_box
      box = [box; box_buffer];
    end
    if is_vel
      vel = [vel; vel_buffer];
    end
    iblock = 1;
  end
end

if iblock > 1
  if is_trj
    trj = [trj; trj_buffer(1:(iblock-1),:)];
  end
  if is_box
    box = [box; box_buffer(1:(iblock-1),:)];
  end
  if is_vel
    vel = [vel; vel_buffer(1:(iblock-1),:)];
  end
end

