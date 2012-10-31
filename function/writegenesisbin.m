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

%% initialization
natom3 = numel(crd);
natom = natom3 / 3;

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


