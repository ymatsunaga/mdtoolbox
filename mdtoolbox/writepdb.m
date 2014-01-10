function writepdb(filename, pdb, index, format_type)
%% writepdb
% write Protein Data Bank (PDB) file
%
%% Syntax
%# writepdb(filename, pdb);
%# writepdb(filename, pdb, index);
%# writepdb(filename, pdb, index, format_type);
%# writepdb(filename, pdb, [], format_type);
%
%% Description
% This code writes only just ATOM or HETATM records. 
% Note that icode fields (code for insetion of residues, see
% "References" below) are ignored in this routine. 
%
% * filename  - input filename of PDB
% * pdb       - structure data
%          record: 'ATOM  ' or 'HETATM' [natomx6 char]
%          serial: Atom serial number. [natomx1 double]
%            name: Atom name. [natomx4 char]
%          altloc: Alternate location indicator. [natomx1 char]
%         resname: Residue name. [natomx3 char]
%         chainid: Chain identifier. [natomx1 char]
%          resseq: Residue sequence number. [natomx1 double]
%           icode: Code for insertion of residues. [natomx1 char] ()
%             xyz: Cartesian coordinate of atom in Angstrom [natomx3 double]
%       occupancy: Occupancy. [natomx1 double]
%      tempfactor: Temperature factor. [natomx1 double]
%         element: Element symbol, right-justified. [natomx2 char]
%          charge: Charge on the atom. [natomx2 char]
% * index       -  index of atoms to be written
% * format_type -  format type [chars. only 'vmd' can be available,
%                  otherwise default(standard?) format is used]
%
%% Example
%# pdb = readpsf('jac.pdb');
%# pdb.xyz(:, 1) = pdb.xyz(:, 1) + 2.5; %translate in x-axis by 2.5 Angstrom
%# writepdb('jac_x.pdb', pdb);
%# writepdb('jac_x_vmd.pdb', pdb, [], 'vmd');
%
%% See also
% readpdb
% 
%% References
% http://www.wwpdb.org/documentation/format33/sect9.html
% ATOM Record Format
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
%  1 -  6        Record name   "ATOM  " or "HETATM"
%  7 - 11        Integer       serial       Atom serial number.
% 13 - 16        Atom          name         Atom name.
% 17             Character     altLoc       Alternate location indicator.
% 18 - 20        Residue name  resName      Residue name.
% 22             Character     chainID      Chain identifier.
% 23 - 26        Integer       resSeq       Residue sequence number.
% 27             AChar         iCode        Code for insertion of residues.
% 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     occupancy    Occupancy.
% 61 - 66        Real(6.2)     tempFactor   Temperature factor.
% 77 - 78        LString(2)    element      Element symbol, right-justified.
% 79 - 80        LString(2)    charge       Charge on the atom.

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% preparation
natom = size(pdb.record, 1);

if (nargin < 3) || (numel(index) == 0)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
end
if iscolumn(index)
  index = index';
end
  
if nargin < 4
  format_type = 'default';
end

%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write file
if strncmpi(format_type, 'vmd', numel('vmd'))
  % VMD format
  for iatom = index
    fprintf(fid, '%6s', pdb.record(iatom, :));
    fprintf(fid, '%5d', mod(pdb.serial(iatom), 100000));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%4s', pdb.name(iatom, :));
    fprintf(fid, '%1s', pdb.altloc(iatom, :));
    %fprintf(fid, '%3s', pdb.resname(iatom, :));
    %fprintf(fid, '%1s', ' ');
    fprintf(fid, '%4s', pdb.resname(iatom, :));
    fprintf(fid, '%1s', pdb.chainid(iatom, :));
    fprintf(fid, '%4d', mod(pdb.resseq(iatom), 10000));
    %fprintf(fid, '%1s', pdb.icode(iatom, :));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
    fprintf(fid, '%6.2f', pdb.occupancy(iatom));
    fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%2s', pdb.element(iatom, :));
    fprintf(fid, '%2s', pdb.charge(iatom, :));
    fprintf(fid, '\n');
  end

elseif strncmpi(format_type, 'namd', numel('namd'))
  % NAMD format
  fprintf(fid, 'CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n');
  for iatom = index
    fprintf(fid, '%6s', pdb.record(iatom, :));
    if pdb.serial(iatom) < 100000
      fprintf(fid, '%5d', pdb.serial(iatom));
    else
      fprintf(fid, '%5s', lower(dec2hex(pdb.serial(iatom))));
    end
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%4s', pdb.name(iatom, :));
    fprintf(fid, '%1s', pdb.altloc(iatom, :));
    %fprintf(fid, '%3s', pdb.resname(iatom, :));
    %fprintf(fid, '%1s', ' ');
    fprintf(fid, '%4s', pdb.resname(iatom, :));
    fprintf(fid, '%1s', pdb.chainid(iatom, :));
    fprintf(fid, '%4d', mod(pdb.resseq(iatom), 10000));
    %fprintf(fid, '%1s', pdb.icode(iatom, :));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
    fprintf(fid, '%6.2f', pdb.occupancy(iatom));
    fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    %fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', pdb.chainid(iatom, :));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%2s', pdb.element(iatom, :));
    %fprintf(fid, '%2s', pdb.charge(iatom, :));
    fprintf(fid, '\n');
  end

else
  % PDB format with some extentions for large digits of atom_serial and resseq
  for iatom = index
    fprintf(fid, '%6s', pdb.record(iatom, :));
    if pdb.serial(iatom) < 100000
      fprintf(fid, '%5d', pdb.serial(iatom));
      fprintf(fid, '%1s', ' ');
    else
      fprintf(fid, '%6d', pdb.serial(iatom));
    end
    fprintf(fid, '%4s', pdb.name(iatom, :));
    fprintf(fid, '%1s', pdb.altloc(iatom, :));
    %fprintf(fid, '%3s', pdb.resname(iatom, :));
    %fprintf(fid, '%1s', ' ');
    fprintf(fid, '%4s', pdb.resname(iatom, :));
    fprintf(fid, '%1s', pdb.chainid(iatom, :));
    if pdb.resseq(iatom) < 10000
      fprintf(fid, '%4d', pdb.resseq(iatom));
      %fprintf(fid, '%1s', pdb.icode(iatom, :));
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
    elseif pdb.resseq(iatom) < 100000
      fprintf(fid, '%5d', pdb.resseq(iatom));
      fprintf(fid, '%1s', ' ');
    else
      fprintf(fid, '%6d', pdb.resseq(iatom));
    end
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
    fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
    fprintf(fid, '%6.2f', pdb.occupancy(iatom));
    fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%1s', ' ');
    fprintf(fid, '%2s', pdb.element(iatom, :));
    fprintf(fid, '%2s', pdb.charge(iatom, :));
    fprintf(fid, '\n');
  end

end

fprintf(fid, 'END\n');


