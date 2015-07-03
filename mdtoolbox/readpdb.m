function [pdb, crd] = readpdb(filename)
%% readpdb
% read Protein Data Bank (PDB) file
%
%% Syntax
%# pdb = readpdb(filename);
%# [pdb, crd] = readpdb(filename);
%# [~, crd] = readpdb(filename);
%
%% Description
% This code only reads ATOM or HETATM records of PDB file.
% Note that icode fields (code for insetion of residues, see
% "References" below) are ignored in this routine. 
%
% * filename  - input filename of PDB
% * pdb       - structure data
%          record: 'ATOM  ' or 'HETATM' [natom x 6 char]
%          serial: Atom serial number. [natom x 1 double]
%            name: Atom name. [natom x 4 char]
%          altloc: Alternate location indicator. [natom x 1 char]
%         resname: Residue name. [natom x 3 char]
%         chainid: Chain identifier. [natom x 1 char]
%          resseq: Residue sequence number. [natom x 1 double]
%           icode: Code for insertion of residues. [natom x 1 char]
%             xyz: Cartesian coordinate of atom in Angstrom [natom x 3 double]
%       occupancy: Occupancy. [natom x 1 double]
%      tempfactor: Temperature factor. [natom x 1 double]
%         element: Element symbol, right-justified. [natom x 2 char]
%          charge: Charge on the atom. [natom x 2 char]
% * crd       - coordinates [1 x natom3 double]
%               if the PDB file contains multiple MODELs
%               coordinates [nmodel x natom3 double]
%         
%% Example
%# pdb = readpsf('jac.pdb');
%
%% See also
% writepdb
% 
%% References
% http://www.wwpdb.org/documentation/format33/sect9.html
% ATOM Record Format
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
%  1 -  6        Record name   "ATOM  " or "HETATM"
%  7 - 11        Integer       serial       Atom serial number. ### we use 7-12 columns ###
% 13 - 16        Atom          name         Atom name.
% 17             Character     altLoc       Alternate location indicator.
% 18 - 20        Residue name  resName      Residue name. ### we use 18-21 columns ###
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
%
% example:
% 123456789012345678901234567890123456789012345678901234567890123456
% ATOM      1  N   MET A   1       9.984  -7.396  -9.702  1.00  1.00
% ATOM  480033 OH2 TIP3W1770     -51.346  17.802  41.770  1.00  1.00
% 

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% read the whole lines into a cell array
%lines = textscan(fid, '%s', 'delimiter', '\n');
lines = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); %for octave
lines = lines{1};
nline = numel(lines);

%% check the number of atoms and multi-models
lines_clean = {};
icount = 0;
nmodel = 0;
model = {};

for iline = 1:nline
  line = lines{iline};
  if strncmpi(line, 'ATOM', numel('ATOM')) || strncmpi(line, 'HETATM', numel('HETATM'))
    icount = icount + 1;
    lines_clean{icount} = line;
  end
  if strncmpi(line, 'ENDMDL', numel('ENDMDL')) || (iline == nline & icount ~= 0)
    natom = icount;
    nmodel = nmodel + 1;
    model{nmodel} = lines_clean;
    icount = 0;
  end
end

%% allocate variables
pdb.record     = repmat('123456', natom, 1);
pdb.serial     = zeros(natom, 1);
pdb.name       = repmat('3456', natom, 1);
pdb.altloc     = repmat('7', natom, 1);
%pdb.resname    = repmat('890', natom, 1);
pdb.resname    = repmat('8901', natom, 1);
pdb.chainid    = repmat('2', natom, 1);
pdb.resseq     = zeros(natom, 1);
%pdb.icode      = repmat('7', natom, 1);
pdb.xyz        = zeros(natom, 3);
pdb.occupancy  = zeros(natom, 1);
pdb.tempfactor = zeros(natom, 1);
pdb.element    = repmat('78', natom, 1);
pdb.charge     = repmat('90', natom, 1);

%% read PDB data
lines = model{1};
for iatom = 1:natom
  iatom;
  line = repmat(' ', 1, 80);
  line(1:length(lines{iatom})) = lines{iatom};
  pdb.record(iatom, :)     = line(1:6);
  %pdb.serial(iatom)        = str2num(line(7:11));
  num = str2num(line(7:12));
  if isempty(num);
    if iatom > 1
      pdb.serial(iatom)    = pdb.serial(iatom-1) + 1;
    else
      pdb.serial(iatom)    = 1;
    end
  else
    pdb.serial(iatom)      = num;
  end
  pdb.name(iatom, :)       = line(13:16);
  pdb.altloc(iatom, :)     = line(17);
  %pdb.resname(iatom, :)    = line(18:20);
  pdb.resname(iatom, :)    = line(18:21);
  pdb.chainid(iatom, :)    = line(22);
  %pdb.resseq(iatom, :)     = str2num(line(23:26));
  %pdb.icode(iatom, :)      = line(27);
  pdb.resseq(iatom, :)     = str2num(line(23:28));
  pdb.xyz(iatom, 1)        = str2num(line(31:38));
  pdb.xyz(iatom, 2)        = str2num(line(39:46));
  pdb.xyz(iatom, 3)        = str2num(line(47:54));
  pdb.occupancy(iatom)     = str2num(line(55:60));
  pdb.tempfactor(iatom)    = str2num(line(61:66));
  pdb.element(iatom, :)    = line(77:78);
  pdb.charge(iatom, :)     = line(79:80);
end

if nargout >= 2
  crd = zeros(nmodel, 3*natom);
  tmp = pdb.xyz';
  crd(1, :) = tmp(:)';
  
  for imodel = 2:nmodel
    xyz = zeros(natom, 3);
    lines = model{imodel};
    for iatom = 1:natom
      line = repmat(' ', 1, 80);
      line(1:length(lines{iatom})) = lines{iatom};
      xyz(iatom, 1) = str2num(line(31:38));
      xyz(iatom, 2) = str2num(line(39:46));
      xyz(iatom, 3) = str2num(line(47:54));
    end
    tmp = xyz';
    crd(imodel, :) = tmp(:)';
  end
end

