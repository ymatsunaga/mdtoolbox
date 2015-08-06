function [trj, box, vel] = readmarbletrjbin(filename, index_atom, index_time)
%% readmarbletrjbin
% read marble binary-format trajectory file
%
%% Syntax
%# trj = readmarbletrjbin(filename);
%# trj = readmarbletrjbin(filename, index_atom);
%# [trj, box] = readmarbletrjbin(filename, index_atom);
%# [trj, box, vel] = readmarbletrjbin(filename, index_atom);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nframe' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input marble trajectory filename
% * index_atom - atom index or logical index specifying atoms to be read
% * trj        - trajectory [nframex3natom double]
% * box        - box size [nframex3 double]
% * vel        - velocities [nframex3natom double]
%
%% Example
%# trj = readmarbletrj('eq.trj');
%
%% See alo
% readmarbletrj
% readmarblecrd
%

%% initialization
if ~exist('index_time', 'var') || isempty(index_time)
  index_time = [];
end

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'r', 'l');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% check the size of file
%rewind
fseek(fid, 0, 'bof');

fseek(fid, 0, 'eof');
eof = ftell(fid);
fseek(fid, 0, 'bof');

%% read header
TRJ_BIN_HEADER = 'MD_SYSTEM_TRJ_OUT';
s           = fread(fid, length(TRJ_BIN_HEADER)+1, 'char');
trj_version = fread(fid, 2, 'int');
trj_flag    = fread(fid, 1, 'int');
natom       = fread(fid, 1, 'int');
natom3      = natom*3;

if ~exist('index_atom', 'var') || isempty(index_atom)
  index_atom = 1:natom;
else
  if islogical(index_atom)
    index_atom = find(index_atom);
  end
end
index_atom3 = to3(index_atom);

%% allocation
trj = zeros(1, numel(index_atom3));
vel = zeros(1, numel(index_atom3));
box = zeros(1, 3);

%% read data
iframe = 0;
idata = 0;
crd  = zeros(1, natom3);
vcrd = zeros(1, natom3);
bcrd = zeros(1, 9);
while 1
  cof = ftell(fid);
  if eof == cof
    break;
  end

  iframe = iframe + 1;
  if isempty(index_time) || ismember(iframe, index_time)
    idata = idata + 1;
    crd  = fread(fid, natom3, 'float64');
    vcrd = fread(fid, natom3, 'float64');
    bcrd = fread(fid, 9, 'float64');
    trj(idata, :) = crd(index_atom3);
    vel(idata, :) = vcrd(index_atom3);
    box(idata, :) = bcrd([1 5 9])';
  else
    fseek(fid, 8*natom3, 0);
    fseek(fid, 8*natom3, 0);
    fseek(fid, 8*9, 0);
  end

  if (~isempty(index_time)) && (iframe >= max(index_time))
    break;
  end
end


