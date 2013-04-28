function writepcrst(filename, pmode, r0, k0, r0a, k0a, nstep1, nstep2)
%% writepcrst
% write configuration file for principal component restraint as a mu2lib input file
%
%% Syntax
%# writepcrst(filename, pmode, r0, r0a, k0, k0a, nstep1, nstep2)
%
%% Description
%
% * filename  - input amber coordinate/restart filename
% * pmode     - coordinates [natom3 x nmode double]
% * r0        - [nmode x 1 double]
% * r0a       - [nmode x 1 double]
% * k0        - [nmode x 1 double]
% * k0a       - [nmode x 1 double]
% * nstep1    - [nmode x 1 double]
% * nstep2    - [nmode x 1 double]
%
%% Example
%
%% See alo
%
%% References
% 

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
nmode = size(pmode, 2);
natom3 = size(pmode, 1);
natom = natom3 / 3;

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write data
fprintf(fid, '# number of modes = %d\n', nmode);
fprintf(fid, '#\n');

if nargin == 4
  for imode = 1:nmode
    fprintf(fid, ' &pcrst\n');
    fprintf(fid, '   r0  = %15.11f, k0  = %13.5f,\n', r0(imode), k0(imode));
    fprintf(fid, '   pcs = ');
    for i = 1:natom3
      fprintf(fid, '%15.11f, ', pmode(i, imode));
      if mod(i, 10) == 0
        fprintf(fid, '\n         ');
      end
    end
    if mod(natom3, 10) ~= 0
      fprintf(fid, '\n');
    end
    fprintf(fid, ' /\n');
  end

else
  for imode = 1:nmode
    fprintf(fid, ' &pcrst\n');
    fprintf(fid, '   ifvari = 1, nstep1 = %d, nstep2 = %d,\n', nstep1(imode), nstep2(imode));
    fprintf(fid, '   r0  = %15.11f, k0  = %13.5f,\n', r0(imode), k0(imode));
    fprintf(fid, '   r0a = %15.11f, k0a = %13.5f,\n', r0a(imode), k0a(imode));
    fprintf(fid, '   pcs = ');
    for i = 1:natom3
      fprintf(fid, '%15.11f, ', pmode(i, imode));
      if mod(i, 10) == 0
        fprintf(fid, '\n         ');
      end
    end
    if mod(natom3, 10) ~= 0
      fprintf(fid, '\n');
    end
    fprintf(fid, ' /\n');
  end
end

