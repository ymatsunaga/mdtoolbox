function writegpr(filename, pdb)
%% writegpr
% write Ca-based Go-model parameter file for GENESIS
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

index = selectname(pdb.name, 'CA') & selectname(pdb.chainid, 'A');
pdb_ca = substruct(pdb, index);
crd_ca = pdb_ca.xyz';
crd_ca = crd_ca(:)';
natom = size(pdb_ca.xyz, 1);

%% calc bond
bond_list = [(1:(natom-1))' (2:natom)'];
bond_value = zeros(size(bond_list, 1), 1);
for i = 1:size(bond_list, 1)
  index = bond_list(i, :);
  bond_value(i) = calcbond(crd_ca(to3(index)));
end

%% calc angle
angle_list = [(1:(natom-2))' (2:(natom-1))' (3:natom)'];
angle_value = zeros(size(angle_list, 1), 1);
for i = 1:size(angle_list, 1)
  index = angle_list(i, :);
  angle_value(i) = calcangle(crd_ca(to3(index))).*180./pi;
end

%% calc dihedral
dihedral_list = [(1:(natom-3))' (2:(natom-2))' (3:(natom-1))' (4:natom)'];
dihedral_value = zeros(size(dihedral_list, 1), 1);
for i = 1:size(dihedral_list, 1)
  index = dihedral_list(i, :);
  dihedral_value(i) = calcdihedral(crd_ca(to3(index)));
end

%% calc native contacts
index = selectname(pdb.record, 'ATOM');
pdb = substruct(pdb, index);
crd = pdb.xyz';
crd = crd(:)';
[contact_list, contact_value] = calccontact(crd, pdb.name, pdb.resseq);

%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% write data
% bond
fprintf(fid, '%8d%s\n', size(bond_list, 1), ' !BOND');
fprintf(fid, '%6d%6d  %13.5e%13.5e\n', ... 
              [bond_list(:,1)'; bond_list(:,2)'; ...
               bond_value'; 100*ones(size(bond_list, 1), 1)']);
fprintf(fid, '\n');

% angle
fprintf(fid, '%8d%s\n', size(angle_list, 1), ' !ANGL');
fprintf(fid, '%6d%6d%6d  %13.5e%13.5e\n', ... 
              [angle_list(:,1)'; angle_list(:,2)'; angle_list(:,3)'; ...
               angle_value'; 20*ones(size(angle_list, 1), 1)']);
fprintf(fid, '\n');

% dihedral
fprintf(fid, '%8d%s\n', size(dihedral_list, 1)*2, ' !DIHE');
for i = 1:size(dihedral_list, 1)
  periodicity = repmat([1 3], 1, size(dihedral_list, 1));
  fprintf(fid, '%6d%6d%6d%6d  %13.5e%13.5e%2d\n', ... 
              [dihedral_list(i,:) (dihedral_value(i).*180./pi + 180.0) 1.0 1]);
  fprintf(fid, '%6d%6d%6d%6d  %13.5e%13.5e%2d\n', ... 
              [dihedral_list(i,:) (dihedral_value(i).*3.*180./pi + 180.0) 0.5 3]);
end
fprintf(fid, '\n');

% improper
fprintf(fid, '%8d%s\n', 0, ' !IMPR');
fprintf(fid, '\n');
fprintf(fid, '\n');

% pair
fprintf(fid, '%8d%s\n', size(contact_list, 1), ' !PAIR');
fprintf(fid, '%6d%6d  %13.5e%13.5e\n', ... 
              [contact_list(:,1)'; contact_list(:,2)'; ...
               ones(size(contact_list, 1), 1)'; contact_value']);
fprintf(fid, '\n');

% non-pair
fprintf(fid, '%8d%s\n', size(dihedral_list, 1), ' !NONPAIR');
fprintf(fid, '%13.5e%13.5e\n', 0.2, 4.0);
fprintf(fid, '\n');


function s = calcdihedral(x)
nframe = size(x,1);
s = zeros(nframe,1);

for iframe = 1:nframe
  d1 = x(iframe,4:6) - x(iframe,1:3);
  d2 = x(iframe,7:9) - x(iframe,4:6);
  d3 = x(iframe,10:12) - x(iframe,7:9);
  m1 = cross(d1,d2);
  m2 = cross(d2,d3);
  s(iframe) = acos(dot(m1,m2)./(norm(m1).*norm(m2)));
  rotdirection = dot(d2,cross(m1,m2));
  if rotdirection < 0
    s(iframe) = 2*pi - s(iframe);
  end
end

