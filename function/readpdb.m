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
%          record: 'ATOM  ' or 'HETATM' [natomx6 char]
%          serial: Atom serial number. [natomx1 double]
%            name: Atom name. [natomx4 char]
%          altloc: Alternate location indicator. [natomx1 char]
%         resname: Residue name. [natomx3 char]
%         chainid: Chain identifier. [natomx1 char]
%          resseq: Residue sequence number. [natomx1 double]
%           icode: Code for insertion of residues. [natomx1 char]
%             xyz: Cartesian coordinate of atom in Angstrom [natomx3 double]
%       occupancy: Occupancy. [natomx1 double]
%      tempfactor: Temperature factor. [natomx1 double]
%         element: Element symbol, right-justified. [natomx2 char]
%          charge: Charge on the atom. [natomx2 char]
% * crd       - coordinates [1 x natom3 double]
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
% e.g.,
% 123456789012345678901234567890123456789012345678901234567890123456
% ATOM      1  N   MET A   1       9.984  -7.396  -9.702  1.00  1.00
% ATOM  480033 OH2 TIP3W1770     -51.346  17.802  41.770  1.00  1.00
% 

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse file
natom = 0;
while ~feof(fid)
  line_raw = strtrim(fgetl(fid));
  
  if strncmpi(line_raw, 'ATOM', numel('ATOM')) | strncmpi(line_raw, 'HETATM', numel('HETATM'))
    % ATOM or HETATM record
    natom = natom + 1;
  end

end
frewind(fid);

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

natom = 0;
while ~feof(fid)
  line_raw = strtrim(fgetl(fid));
  
  if strncmpi(line_raw, 'ATOM', numel('ATOM')) | strncmpi(line_raw, 'HETATM', numel('HETATM'))
    % ATOM or HETATM record
    line = repmat(' ', 1, 80);
    line(1:length(line_raw)) = line_raw;
    natom = natom + 1;
    pdb.record(natom, :)     = line(1:6);
    %pdb.serial(natom)        = str2num(line(7:11));
    num = str2num(line(7:12));
    if numel(num) == 0;
      if natom > 1
        pdb.serial(natom)    = pdb.serial(natom-1) + 1;
      else
        pdb.serial(natom)    = 1;
      end
    else
      pdb.serial(natom)      = num;
    end
    pdb.name(natom, :)       = line(13:16);
    pdb.altloc(natom, :)     = line(17);
    %pdb.resname(natom, :)    = line(18:20);
    pdb.resname(natom, :)    = line(18:21);
    pdb.chainid(natom, :)    = line(22);
    %pdb.resseq(natom, :)     = str2num(line(23:26));
    %pdb.icode(natom, :)      = line(27);
    pdb.resseq(natom, :)     = str2num(line(23:28));
    pdb.xyz(natom, 1)        = str2num(line(31:38));
    pdb.xyz(natom, 2)        = str2num(line(39:46));
    pdb.xyz(natom, 3)        = str2num(line(47:54));
    pdb.occupancy(natom)     = str2num(line(55:60));
    pdb.tempfactor(natom)    = str2num(line(61:66));
    pdb.element(natom, :)    = line(77:78);
    pdb.charge(natom, :)     = line(79:80);
  end

end

if nargout >= 2
  crd = pdb.xyz';
  crd = crd(:)';
end

