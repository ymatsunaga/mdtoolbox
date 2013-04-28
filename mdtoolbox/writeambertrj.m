function writeambertrj(filename, trj, box, title)
%% writeambertrj
% write amber ascii-format trajectory format file
%
%% Syntax
%# writeambertrj(filename, trj)
%# writeambertrj(filename, trj, box)
%# writeambertrj(filename, trj, box, title)
%# writeambertrj(filename, trj, [], title)
%
%% Description
% This code puts trajectories into 
% an amber trajectory format file. 
% If box information is given, 
% box sizes are appended.
%
% * filename  - output dcd trajectory filename
% * trj       - trajectory [nstep x natom3 double]
% * box       - box size [nstep x 3 double]
% * title     - title characters [chars]
%
%% Example
%# natom = 3343;
%# trj = readambertrj(natom, '4ake.trj');
%# trj(:, 1:3:end) = trj(:, 1:3:end) + 1.5;
%# writeambertrj('4ake_translated.trj', trj, [], 'translated in x axis')
%
%% See also
% readambertrj
% readambertrjbox
%
%% References
% http://ambermd.org/formats.html#trajectory
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
natom3 = size(trj, 2);
nstep = size(trj, 1);

if nargin < 4 | (isempty(title))
  title = sprintf('FILENAME=%s CREATED BY MATLAB', filename);
end

%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write data
fprintf(fid, '%s\n', title);
for istep = 1:nstep
  for i = 1:10:natom3
    fprintf(fid, '%8.3f', trj(istep, i:min(i+9,natom3)));
    fprintf(fid, '\n');
  end
  if (nargin >= 3) & (~isempty(box))
    fprintf(fid, '%8.3f', box(istep, :));
    fprintf(fid, '\n');
  end
end

