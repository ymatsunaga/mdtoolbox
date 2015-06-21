function parm = readprmtop(filename)
%% readprmtop
% read amber parameter/topology file
%
%% Syntax
%# prmtop = readprmtop(filename)
%
%% Description
% Output is a strucure variable. 
%
% * fname     - input filename of prmtop
% * prmtop    - structure data 
%                         title: [20x4 char]
%                         natom: double
%                        ntypes: double
%                         nbonh: double
%                         mbona: double
%                        ntheth: double
%                        mtheta: double
%                         nphih: double
%                         mphia: double
%                        nhparm: double
%                         nparm: double
%                           nnb: double
%                          nres: double
%                         nbona: double
%                        ntheta: double
%                         nphia: double
%                        numbnd: double
%                        numang: double
%                         nptra: double
%                         natyp: double
%                          nphb: double
%                        ifpert: double
%                         nbper: double
%                         ngper: double
%                         ndper: double
%                         mbper: double
%                         mgper: double
%                         mdper: double
%                         ifbox: double
%                         nmxrs: double
%                         ifcap: double
%                      numextra: double
%                     atom_name: [natomx4 char]
%                        charge: [natomx1 double]
%                          mass: [natomx1 double]
%               atom_type_index: [natomx1 double]
%         number_excluded_atoms: [natomx1 double]
%          nonbonded_parm_index: [XXXx1 double]
%                 residue_label: [nresx4 char]
%               residue_pointer: [nresx1 double]
%           bond_force_constant: [numbndx1 double]
%              bond_equil_value: [numbndx1 double]
%          angle_force_constant: [numangx1 double]
%             angle_equil_value: [numangx1 double]
%       dihedral_force_constant: [nptrax1 double]
%          dihedral_periodicity: [nptrax1 double]
%                dihedral_phase: [nptrax1 double]
%                         solty: [natypx1 double]
%           lennard_jones_acoef: [nlennx1 double]
%           lennard_jones_bcoef: [nlennx1 double]
%            bonds_inc_hydrogen: [nhydrox3 double]
%        bonds_without_hydrogen: [nonhydrox3 double]
%           angles_inc_hydrogen: [nanglehydrox4 double]
%       angles_without_hydrogen: [nanglenonhydrox4 double]
%        dihedrals_inc_hydrogen: [ndihedralhydrox5 double]
%    dihedrals_without_hydrogen: [ndihedralnonhydrox5 double]
%           excluded_atoms_list: [nexcludex1 double]
%                   hbond_acoef: []
%                   hbond_bcoef: []
%                         hbcut: []
%               amber_atom_type: [natomx4 char]
%     tree_chain_classification: [natomx4 char]
%                    join_array: [natomx1 double]
%                        irotat: [natomx1 double]
%                    radius_set: [1x80 char]
%                         radii: [natomx1 double]
%                        screen: [natomx1 double]
% 
%% Example
%# prmtop = readpsf('ak.prmtop');
%
%% See also
% readpsf
% 

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse file
parm = struct;

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if strncmpi(line, '%FLAG ', numel('%FLAG '))
    key = line((numel('%FLAG ')+1):end);
    key = strtrim(key);
    key = lower(key);

    while ~strncmpi(line, '%FORMAT', numel('%FORMAT'))
      line = strtrim(fgetl(fid));
    end
    s_format = regexp(line, '\(.*\)', 'match');
    s_format = s_format{1}(2:(end-1));

    data = [];

    if regexpi(s_format, 'E')
      s_format = regexpi(s_format, 'E', 'split');
      s_format = s_format{end};
      if regexp(s_format, '\.')
        s_format = regexpi(s_format, '\.', 'split');
        s_format = s_format{1};
      end
      s_format = ['%' s_format 'e'];
      data = fscanf(fid, s_format);
    
    elseif regexpi(s_format, 'F')
      s_format = regexpi(s_format, 'F', 'split');
      s_format = s_format{end};
      if regexp(s_format, '\.')
        s_format = regexpi(s_format, '\.', 'split');
        s_format = s_format{1};
      end
      s_format = ['%' s_format 'f'];
      data = fscanf(fid, s_format);
    
    elseif regexpi(s_format, 'I')
      s_format = regexpi(s_format, 'I', 'split');
      s_format = s_format{end};
      s_format = ['%' s_format 'd'];
      data = fscanf(fid, s_format);
    
    elseif regexpi(s_format, 'a')
      s_format = regexpi(s_format, 'a', 'split');
      s_format = s_format{end};
      data = fscanf(fid, '%[^%]\n');
      newline = sprintf('\n');
      data = strrep(data, newline, '');
      data = reshape(data, str2num(s_format), numel(data)/str2num(s_format));
      data = data';
      
    else
      disp(line);
      disp(key);
      disp(s_format);
      error('could not understand the format of flag');
    
    end
    
    if strncmpi(key, 'BONDS_INC_HYDROGEN', numel('BONDS_INC_HYDROGEN'))
      data = reshape(data, 3, [])';
    elseif strncmpi(key, 'BONDS_WITHOUT_HYDROGEN', numel('BONDS_WITHOUT_HYDROGEN'))
      data = reshape(data, 3, [])';
    elseif strncmpi(key, 'ANGLES_INC_HYDROGEN', numel('ANGLES_INC_HYDROGEN'))
      data = reshape(data, 4, [])';
    elseif strncmpi(key, 'ANGLES_WITHOUT_HYDROGEN', numel('ANGLES_WITHOUT_HYDROGEN'))
      data = reshape(data, 4, [])';
    elseif strncmpi(key, 'DIHEDRALS_INC_HYDROGEN', numel('DIHEDRALS_INC_HYDROGEN'))
      data = reshape(data, 5, [])';
    elseif strncmpi(key, 'DIHEDRALS_WITHOUT_HYDROGEN', numel('DIHEDRALS_WITHOUT_HYDROGEN'))
      data = reshape(data, 5, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_PARAMETER_01', numel('CHARMM_CMAP_PARAMETER_01'))
      data = reshape(data, 8, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_PARAMETER_02', numel('CHARMM_CMAP_PARAMETER_02'))
      data = reshape(data, 8, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_PARAMETER_03', numel('CHARMM_CMAP_PARAMETER_03'))
      data = reshape(data, 8, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_PARAMETER_04', numel('CHARMM_CMAP_PARAMETER_04'))
      data = reshape(data, 8, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_PARAMETER_05', numel('CHARMM_CMAP_PARAMETER_05'))
      data = reshape(data, 8, [])';
    elseif strncmpi(key, 'CHARMM_CMAP_INDEX', numel('CHARMM_CMAP_INDEX'))
      data = reshape(data, 6, [])';
    end
      
    if strncmpi(key, 'POINTERS', numel('POINTERS'))
      i = 1;
      if i <= numel(data); parm.natom = data(i); end; i=i+1;    % total number of atoms 
      if i <= numel(data); parm.ntypes = data(i); end; i=i+1;   % total number of distinct atom types
      if i <= numel(data); parm.nbonh = data(i); end; i=i+1;    % number of bonds containing hydrogen
      if i <= numel(data); parm.mbona = data(i); end; i=i+1;    % number of bonds not containing hydrogen
      if i <= numel(data); parm.ntheth = data(i); end; i=i+1;   % number of angles containing hydrogen
      if i <= numel(data); parm.mtheta = data(i); end; i=i+1;   % number of angles not containing hydrogen
      if i <= numel(data); parm.nphih = data(i); end; i=i+1;    % number of dihedrals containing hydrogen
      if i <= numel(data); parm.mphia = data(i); end; i=i+1;    % number of dihedrals not containing hydrogen
      if i <= numel(data); parm.nhparm = data(i); end; i=i+1;   % currently not used
      if i <= numel(data); parm.nparm = data(i); end; i=i+1;    % used to determine if addles created prmtop
      if i <= numel(data); parm.nnb = data(i); end; i=i+1;      % number of excluded atoms
      if i <= numel(data); parm.nres = data(i); end; i=i+1;     % number of residues
      if i <= numel(data); parm.nbona = data(i); end; i=i+1;    % mbona + number of constraint bonds
      if i <= numel(data); parm.ntheta = data(i); end; i=i+1;   % mtheta + number of constraint angles
      if i <= numel(data); parm.nphia = data(i); end; i=i+1;    % mphia + number of constraint dihedrals
      if i <= numel(data); parm.numbnd = data(i); end; i=i+1;   % number of unique bond types
      if i <= numel(data); parm.numang = data(i); end; i=i+1;   % number of unique angle types
      if i <= numel(data); parm.nptra = data(i); end; i=i+1;    % number of unique dihedral types
      if i <= numel(data); parm.natyp = data(i); end; i=i+1;    % number of atom types in parameter file, see solty below
      if i <= numel(data); parm.nphb = data(i); end; i=i+1;     % number of distinct 10-12 hydrogen bond pair types
      if i <= numel(data); parm.ifpert = data(i); end; i=i+1;   % set to 1 if perturbation info is to be read in
      if i <= numel(data); parm.nbper = data(i); end; i=i+1;    % number of bonds to be perturbed
      if i <= numel(data); parm.ngper = data(i); end; i=i+1;    % number of angles to be perturbed
      if i <= numel(data); parm.ndper = data(i); end; i=i+1;    % number of dihedrals to be perturbed
      if i <= numel(data); parm.mbper = data(i); end; i=i+1;    % number of bonds with atoms completely in perturbed group
      if i <= numel(data); parm.mgper = data(i); end; i=i+1;    % number of angles with atoms completely in perturbed group
      if i <= numel(data); parm.mdper = data(i); end; i=i+1;    % number of dihedrals with atoms completely in perturbed groups
      if i <= numel(data); parm.ifbox = data(i); end; i=i+1;    % set to 1 if standard periodic box, 2 when truncated octahedral
      if i <= numel(data); parm.nmxrs = data(i); end; i=i+1;    % number of atoms in the largest residue
      if i <= numel(data); parm.ifcap = data(i); end; i=i+1;    % set to 1 if the cap option from edit was specified
      if i <= numel(data); parm.numextra = data(i); end; i=i+1; % number of extra points found in topology
      if i <= numel(data); parm.ncopy = data(i); end;           % number of pimd slices / number of beads
    else
      parm = setfield(parm, key, data);
    end
  
  end
  
end


%% convert excluded_atoms_list to a pair list (just for convenience)
parm.excluded_pair = zeros(numel(parm.excluded_atoms_list), 2);
istart = 0;
for iatom = 1:parm.natom
  num = parm.number_excluded_atoms(iatom);
  index = (istart+1):(istart+num);
  parm.excluded_pair(index, :) = [iatom*ones(num, 1) parm.excluded_atoms_list(index)];
  istart = istart + num;
end

index = find(parm.excluded_pair(:,2) == 0);
parm.excluded_pair(index, :) = [];


%% create other elements to be consistent with psf and convenient for atom selection
natom = numel(parm.charge);

% atom_id
parm.atom_id = (1:natom)';

% segment_name
% todo...

% residue_id
parm.residue_id = zeros(natom, 1);
iresidue = 1;
for i = 1:(numel(parm.residue_pointer)-1)
  index = parm.residue_pointer(i):(parm.residue_pointer(i+1) - 1);
  parm.residue_id(index) = iresidue;
  iresidue = iresidue + 1;
end
index = parm.residue_pointer(end):natom;
parm.residue_id(index) = iresidue;

% residue_name
parm.residue_name = repmat('1234', natom, 1);
iresidue = 1;
for i = 1:(numel(parm.residue_pointer)-1)
  index = parm.residue_pointer(i):(parm.residue_pointer(i+1) - 1);
  parm.residue_name(index, :) = repmat(parm.residue_label(i, :), numel(index), 1);
  iresidue = iresidue + 1;
end
index = parm.residue_pointer(end):natom;
parm.residue_name(index, :) = repmat(parm.residue_label(end, :), numel(index), 1);

