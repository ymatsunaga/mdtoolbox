function writepsf(filename, psf)
%% writepsf
% writes charmm or xplor type Protein Structure File (PSF)
%
%% Syntax
%# writepsf(filename, psf)
%
%% Description
%
% * filename  - output filename of psf
% * psf       - structure data 
%           isPSF: logical
%           isEXT: logical
%          isCMAP: logical
%          isCHEQ: logical
%           title: [ntitlex80 char]
%         atom_id: [natomx1 double]
%    segment_name: [natomx4 char]
%      residue_id: [natomx1 double]
%    residue_name: [natomx4 char]
%       atom_name: [natomx4 char]
%       atom_type: [natomx4 char] or [natomx1 integer]
%          charge: [natomx1 double]
%            mass: [natomx1 double]
%           nbond: [double]
%       bond_list: [nbondx2 double]
%          ntheta: [double]
%      theta_list: [nthetax3 double]
%            nphi: [double]
%        phi_list: [nphix4 double]
%          nimphi: [double]
%      imphi_list: [nimphix4 double]
%            ndon: [double]
%        don_list: [ndonx2 double]
%            nacc: [double]
%        acc_list: [naccx2 double]
%             nnb: [double]
%         nb_list: [nnbx1 double]
%            ngrp: [double]
%        grp_list: [ngrpx3 double]
%         ncrterm: [double]
%     crterm_list: [ncrtermx8 double]
% 
%% Example
%# 
%
%% References
% NAMD Tutorial
% http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node21.html
% User Manual for EGO_VIII 
% http://heller.userweb.mwn.de/ego/manual/node88.html
% 
%% See also
% readpsf readcharmmparam
% 

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% header
line = [];
if isfield(psf, 'isPSF')
  if psf.isPSF
    line = [line 'PSF '];
  end
end
if isfield(psf, 'isEXT')
  if psf.isEXT
    line = [line 'EXT '];
  end
else
  psf.isEXT = false;
end
if isfield(psf, 'isCMAP')
  if psf.isCMAP
    line = [line 'CMAP'];
  end
end
if isfield(psf, 'isCHEQ')
  if psf.isCHEQ
    line = [line 'CHEQ'];
  end
end
for i = (numel(line)+1):80
  line = [line ' '];
end
fprintf(fid, '%s\n', line);
fprintf(fid, '\n');

%% format
if psf.isEXT
  fmt_list = '%10d';
  if ischar(psf.atom_type)
    fmt_atom = '%10d %-8s %-8d %8s %8s %4s %14f%14f%12d\n';
  else
    fmt_atom = '%10d %-8s %-8d %8s %8s %4d %14f%14f%12d\n';
  end
else
  fmt_list = '%8d';
  if ischar(psf.atom_type)
    fmt_atom = '%8d %-4s %-4d %4s %4s %4s %14f%14f%12d\n';
  else
    fmt_atom = '%8d %-4s %-4d %4s %4s %4d %14f%14f%12d\n';
  end
end

%% title
ntitle = size(psf.title, 1);
fprintf(fid, [fmt_list ' !NTITLE\n'], ntitle);
for i = 1:ntitle
  fprintf(fid, '%s\n', psf.title(i, :));
end
fprintf(fid, '\n');

%% atom information
natom = numel(psf.atom_id);
fprintf(fid, [fmt_list ' !NATOM\n'], natom);
for i = 1:natom
  fprintf(fid, fmt_atom, psf.atom_id(i), psf.segment_name(i, :), psf.residue_id(i, :), psf.residue_name(i, :), psf.atom_name(i, :), psf.atom_type(i, :), psf.charge(i), psf.mass(i), 0);
end
fprintf(fid, '\n');

%% bonds
fprintf(fid, [fmt_list ' !NBOND: bonds\n'], psf.nbond);
list = psf.bond_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% angles
fprintf(fid, [fmt_list ' !NTHETA: angles\n'], psf.ntheta);
list = psf.theta_list';
n = numel(list);
for i = 1:9:n
  fprintf(fid, fmt_list, list(i:min(i+8,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% dihedrals
fprintf(fid, [fmt_list ' !NPHI: dihedrals\n'], psf.nphi);
list = psf.phi_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% impropers
fprintf(fid, [fmt_list ' !NIMPHI: impropers\n'], psf.nimphi);
list = psf.imphi_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% donors
fprintf(fid, [fmt_list ' !NDON: donors\n'], psf.ndon);
list = psf.don_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% donors
fprintf(fid, [fmt_list ' !NACC: acceptors\n'], psf.nacc);
list = psf.acc_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

%% NNB
%     elseif strncmpi(key(2:end), 'NNB', numel('NNB'))
%       psf.nnb = num;
%       psf.nb_list = fscanf(fid, fmt_list, [1, num]);
%       psf.nb_list = psf.nb_list';
fprintf(fid, [fmt_list ' !NNB\n'], 0);
list = zeros(natom, 1);
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');
    
%% NGRP
%     elseif strncmpi(key(2:end), 'NGRP', numel('NGRP'))
%       psf.ngrp = num;
%       psf.grp_list = fscanf(fid, fmt_list, [3, num]);
%       psf.grp_list = psf.grp_list';

%% cross-terms
fprintf(fid, [fmt_list ' !NCRTERM: cross-terms\n'], psf.ncrterm);
list = psf.crterm_list';
n = numel(list);
for i = 1:8:n
  fprintf(fid, fmt_list, list(i:min(i+7,n)));
  fprintf(fid, '\n');
end
fprintf(fid, '\n');

