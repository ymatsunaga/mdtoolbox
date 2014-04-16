function writegenesisbin(filename, crd, vel, box, header)
%% writegenesisbin
% write genesis binary restart file
%
%% Syntax
%# writegenesisbin(filename, crd, vel, box);
%
%% Description
% write genesis binary restart file
%
% * filename  - filename of genesis biary restart
% * crd       - coordinates [1 x natom3]
% * vel       - velocities [1 x natom3]
% * box       - box size [1 x 3]
%
%% Example
%# writegenesisbin('run.rst', crd, vel, box);
%
%% See also
% readgenesisbin
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
natom3 = numel(crd);
natom = natom3 / 3;

if ~exist('header', 'var') || isempty(header)
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
  header.title1 = title1;
  header.title2 = title2;
  header.nAtom = natom;
  header.rstfile_type = 2;
  header.iseed = 777;
  header.num_deg_freedom = natom*3 - 6;
  header.thermostat_friction = 0;
  header.barostat_friction = [0 0 0];
end

fid = fopen(filename, 'w', 'l');

fwrite(fid, 80, 'int32');
fwrite(fid, header.title1(1:80), '*char');
fwrite(fid, 80, 'int32');
%88

fwrite(fid, 80, 'int32');
fwrite(fid, header.title2(1:80), '*char');
fwrite(fid, 80, 'int32');
%88 + 88

fwrite(fid, 4, 'int32');
fwrite(fid, header.nAtom, 'int32');
fwrite(fid, 4, 'int32');
%88 + 88 + 12

fwrite(fid, 4, 'int32');
fwrite(fid, header.rstfile_type, 'int32');
fwrite(fid, 4, 'int32');
%88 + 88 + 12 + 12

fwrite(fid, 8, 'int32');
fwrite(fid, header.iseed, 'int32');
fwrite(fid, header.num_deg_freedom, 'int32');
fwrite(fid, 8, 'int32');
%88 + 88 + 12 + 12 + 16

fwrite(fid, 24, 'int32');
fwrite(fid, box(1), 'double');
fwrite(fid, box(2), 'double');
fwrite(fid, box(3), 'double');
fwrite(fid, 24, 'int32');
%88 + 88 + 12 + 12 + 16 + 28

fwrite(fid, 8, 'int32');
fwrite(fid, header.thermostat_friction, 'double');
fwrite(fid, 8, 'int32');
%88 + 88 + 12 + 12 + 16 + 28 + 12

fwrite(fid, 24, 'int32');
fwrite(fid, header.barostat_friction, 'double');
fwrite(fid, 24, 'int32');
%88 + 88 + 12 + 12 + 16 + 28 + 12 + 12

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, crd(1:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, crd(2:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, crd(3:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, vel(1:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, vel(2:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fwrite(fid, 8*header.nAtom, 'int32');
fwrite(fid, vel(3:3:end), 'double');
fwrite(fid, 8*header.nAtom, 'int32');

fclose(fid);


