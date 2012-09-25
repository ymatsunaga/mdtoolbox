function psf = readpsf(filename)
%% readpsf
% read charmm or xplor type Protein Structure File (PSF)
%
%% Syntax
%# psf = readpsf(filename);
%
%% Description
% This code supports both charmm-type and xplor-type psf.
% Output is a strucure variable. 
%
% * filename  - input filename of psf
% * psf       - structure data 
%           title: [ntitlex80 char]
%         atom_id: [natomx1 double]
%    segment_name: [natomx4 char]
%      residue_id: [natomx1 double]
%    residue_name: [natomx4 char]
%       atom_name: [natomx4 char]
%       atom_type: [natomx4 char]
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
%# psf = readpsf('jac.psf');
%# psf = readpsf('jac_xplor.psf');
%
%% See also
% readamberparm
% 

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse file
isPSF = false;
isEXT = false;
isCMAP = false;
isCHEQ = false;

line = fgetl(fid);
line_size = numel(line);
for i=1:4:(line_size-3)
  if strcmpi(line(i:(i+3)), 'PSF '); isPSF = true; end;
  if strcmpi(line(i:(i+3)), 'EXT '); isEXT = true; end;
  if strcmpi(line(i:(i+3)), 'CMAP'); isCMAP = true; end;
  if strcmpi(line(i:(i+3)), 'CHEQ'); isCHEQ = true; end;
end

if ~isPSF
  error('Error: sorry, it seems not be a PSF file')
end

if isEXT
  fmt_atom = '%10d %8s %8d %8s %8s %4s %14f%14f%8d%*[^\n]';
  fmt_list = '%10d';
else
  fmt_atom = '%8d %4s %4d %4s %4s %4s %14f%14f%8d%*[^\n]';
  fmt_list = '%8d';
end

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if regexp(line, '.*\d.*!\w.*')
    [num, key] = strtok(line, '!');
    num = str2num(num);
    num = num(1);

    if strncmpi(key(2:end), 'NTITLE', numel('NTITLE'))
      psf.title = '';
      for i = 1:num
        psf.title = strvcat(psf.title, fgetl(fid));
      end
    
    elseif strncmpi(key(2:end), 'NATOM', numel('NATOM'))
      C = textscan(fid, fmt_atom, num);
      psf.atom_id      = double(C{1});
      psf.segment_name = char(C{2});
      psf.residue_id   = double(C{3});
      psf.residue_name = char(C{4});
      psf.atom_name    = char(C{5});
      if isletter(C{6}{1}(1))
        psf.atom_type  = char(C{6});
      else
        psf.atom_type  = str2num(char(C{6}));
      end
      psf.charge       = double(C{7});
      psf.mass         = double(C{8});
      
    elseif strncmpi(key(2:end), 'NBOND', numel('NBOND'))
      psf.nbond = num;
      psf.bond_list = fscanf(fid, fmt_list, [2, num]);
      psf.bond_list = psf.bond_list';
    
    elseif strncmpi(key(2:end), 'NTHETA', numel('NTHETA'))
      psf.ntheta = num;
      psf.theta_list = fscanf(fid, fmt_list, [3, num]);
      psf.theta_list = psf.theta_list';
    
    elseif strncmpi(key(2:end), 'NPHI', numel('NPHI'))
      psf.nphi = num;
      psf.phi_list = fscanf(fid, fmt_list, [4, num]);
      psf.phi_list = psf.phi_list';
    
    elseif strncmpi(key(2:end), 'NIMPHI', numel('NIMPHI'))
      psf.nimphi = num;
      psf.imphi_list = fscanf(fid, fmt_list, [4, num]);
      psf.imphi_list = psf.imphi_list';
    
    elseif strncmpi(key(2:end), 'NDON', numel('NDON'))
      psf.ndon = num;
      psf.don_list = fscanf(fid, fmt_list, [2, num]);
      psf.don_list = psf.don_list';
    
    elseif strncmpi(key(2:end), 'NACC', numel('NACC'))
      psf.nacc = num;
      psf.acc_list = fscanf(fid, fmt_list, [2, num]);
      psf.acc_list = psf.acc_list';
    
    elseif strncmpi(key(2:end), 'NNB', numel('NNB'))
      psf.nnb = num;
      psf.nb_list = fscanf(fid, fmt_list, [1, num]);
      psf.nb_list = psf.nb_list';
    
    elseif strncmpi(key(2:end), 'NGRP', numel('NGRP'))
      psf.ngrp = num;
      psf.grp_list = fscanf(fid, fmt_list, [3, num]);
      psf.grp_list = psf.grp_list';
    
    elseif strncmpi(key(2:end), 'NCRTERM', numel('NCRTERM'))
      psf.ncrterm = num;
      psf.crterm_list = fscanf(fid, fmt_list, [8, num]);
      psf.crterm_list = psf.crterm_list';
    
    end
  end
end


