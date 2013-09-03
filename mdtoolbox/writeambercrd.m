function writeambercrd(filename, crd, box, vel, title)
%% writeambercrd
% write amber coordinate/restart file
%
%% Syntax
%# writeambercrd(filename, crd);
%# writeambercrd(filename, crd, box);
%# writeambercrd(filename, crd, box, vel);
%# writeambercrd(filename, crd, [], vel);
%# writeambercrd(filename, crd, box, vel, title);
%# writeambercrd(filename, crd, [], [], title);
%
%% Description
% This code puts coordinates, velocities,
% and box sizes to an amber coordinates/restart file.
% The variables should have the coordinates(velocities)
% in order [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename  - input amber coordinate/restart filename
% * crd       - coordinates [1 x natom3 double]
% * vel       - velocities [1 x natom3 double]
% * box       - size of the periodic box [nstep x 3 double]
% * title    - title characters [chars]
%
%% Example
%# crd = readambercrd('ak.crd');
%# crd(1, 1:3:end) = crd(1, 1:3:end) + 10;
%# writeambercrd('ak_translated.crd', crd);
%
%% See alo
% readambercrd
%
%% References
% http://ambermd.org/formats.html#restart
% 

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
natom3 = size(crd, 2);
natom = natom3 / 3;

if (nargin < 5) || (isempty(title))
  title = sprintf('FILENAME=%s CREATED BY MATLAB', filename);
end
for i = (numel(title)+1):80
  title = [title ' '];
end

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write data
% title
fprintf(fid, '%s\n', title);

% natom
%fprintf(fid,'%5d\n',natom);
fprintf(fid, '%6d%15.7e\n', natom, 0.0);

% coordinates
for i = 1:6:natom3
  fprintf(fid, '%12.7f', crd(1, i:min(i+5,natom3)));
  fprintf(fid, '\n');
end

% velocities
if (nargin >= 4) && (~isempty(vel))
  for i = 1:6:natom3
    fprintf(fid, '%12.7f', vel(1, i:min(i+5,natom3)));
    fprintf(fid, '\n');
  end
end

% box sizes
if (nargin >= 3) && (~isempty(box))
  fprintf(fid, '%12.7f', box);
  fprintf(fid, '%12.7f', [90.0 90.0 90.0]);
  fprintf(fid, '\n');
end

