function param = readcharmmparam(filename)
%% readpsf
% read charmm parameter file
%
%% Syntax
%# param = readcharmmparam(filename)
%
%% Description
% 
%% Example
%# param = readcharmmparam('charmm.param');
%
%% See also
% readpsf
% 

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse file
% title
param.title = [];
while ~feof(fid)
  line = fgetl(fid);
  if (numel(line) > 0) && (strncmpi(line(1), '*', numel('*')))
    for i = (numel(line)+1):80
      line = [line ' '];
    end
    param.title = [param.title; line(1:80)];
  else
    break;
  end
end

% bond
fseek(fid, 0, 'bof');
while ~feof(fid)
  line = fgetl(fid);
  line = uncomment(line);
  if strncmpi(line, 'BOND', numel('BOND'))
    C = textscan(fid, '%s %s %f %f', 'CommentStyle', '!', 'CollectOutput', true);
    break;
  end
end
index = calcindex(C);
param.bond_type = C{1}(index, :);
param.bond_value = C{2}(index, :);

% angle
fseek(fid, 0, 'bof');
while ~feof(fid)
  line = fgetl(fid);
  line = uncomment(line);
  if strncmpi(line, 'ANGLe', numel('ANGLe'))
    C = textscan(fid, '%s %s %s %f %f', 'CommentStyle', '!', 'CollectOutput', true);
    break;
  end
end
index = calcindex(C);
param.angle_type = C{1}(index, :);
param.angle_value = C{2}(index, :);

% dihedral
fseek(fid, 0, 'bof');
while ~feof(fid)
  line = fgetl(fid);
  line = uncomment(line);
  if strncmpi(line, 'DIHEdral', numel('DIHEdral')) || strncmpi(line, 'PHI', numel('PHI'))
    C = textscan(fid, '%s %s %s %s %f %f %f', 'CommentStyle', '!', 'CollectOutput', true);
    break;
  end
end
index = calcindex(C);
param.dihedral_type = C{1}(index, :);
param.dihedral_value = C{2}(index, :);

% nonbonded
fseek(fid, 0, 'bof');
while ~feof(fid)
  line = fgetl(fid);
  line = uncomment(line);
  if strncmpi(line, 'NONBONDed', numel('NBONDed'))
    break;
  end
end
while ~feof(fid)
  if strncmpi(line(end), '-', numel('-'))
    line = fgetl(fid);
  else
    break
  end
end
C = textscan(fid, '%s %f %f %f', 'CommentStyle', '!', 'CollectOutput', true);
index = calcindex(C);
param.nonbonded_type = C{1}(index, :);
param.nonbonded_value = C{2}(index, :);

% nbfix
fseek(fid, 0, 'bof');
while ~feof(fid)
  line = fgetl(fid);
  line = uncomment(line);
  if strncmpi(line, 'NBFIX', numel('NBFIX'))
    C = textscan(fid, '%s %s %f %f', 'CommentStyle', '!', 'CollectOutput', true);
    break;
  end
end
index = calcindex(C);
param.nbfix_type = C{1}(index, :);
param.nbfix_value = C{2}(index, :);


%% uncomment line
function result = uncomment(s)
result = regexprep(s, '!.*', '', 'once');
result = strtrim(result);

%% calc index of the line which containts extra keys
function index = calcindex(c)
keys = {'BOND', 'ANGLe', 'THETa', 'DIHEdral' 'PHI', 'IMPH', 'CMAP', 'NONBONDed', 'NBFIX', 'HBONDs', 'END'};
type = c{1};
index = 1:size(type,1);
for i = size(type,1):-1:1
  for j = 1:numel(keys)
    if strncmpi(type{i,1}, keys{j}, numel(keys{j}))
      index = 1:(i-1);
    end
  end
end


