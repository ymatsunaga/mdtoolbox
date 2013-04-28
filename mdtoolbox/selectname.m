function logical_index = selectname(namelist, varargin)
%% selectname
% used for atom selection. Finds all the atoms or residues which matches given names.
%
%% Syntax
%# logical_index = selectname(namelist, name)
%# logical_index = selectname(namelist, name, name)
%# logical_index = selectname(namelist, name, name, ...)
%
%% Description
% Used for atom selection. 
% Compares namelist(an array of strings) and input names and 
% returns a logical index which has 'trues' for the matched atoms. 
% Wild cards '*' are available for input name. 
% Multiple input names can be specified. If multiple names are given, 
% boolean 'OR' operations are applied for multiple name matching. 
% Name are compared in case-insensitive manner.
% Blanks(whitespaces) in the front and rear of namelist and input
% names are ignored. 
%
% * namelist - atom name list [natom x n char]
% * name     - selection name [1 x m char]
%
%% Example
%# % read amber parm file
%# parm = readparm('ak.parm');
%#
%# % choose CA atoms
%# logical_index = selectname(parm.atom_name, 'CA');
%# index_ca = find(logical_index);
%#
%# % choose backbone atoms
%# logical_index = selectname(parm.atom_name, 'CA', 'C', 'N', 'O');
%# index_bb = find(logical_index);
%#
%# % choose non-hydrogen atoms
%# logical_index = ~selectname(parm.atom_name, 'H*');
%# index_heavy = find(logical_index);
%

[natom, ~] = size(namelist);
logical_index = false(natom, 1);

for i = 1:numel(varargin)
  %logical_index = logical_index | all(bsxfun(@eq, namelist, varargin{i}), 2);
  query = strtrim(varargin{i});
  if regexp(query, '\*')
    % matching with wild card
    query = ['^' query];
    for iatom = 1:natom
      if regexp(strtrim(namelist(iatom, :)), query)
        logical_index(iatom) = true;
      end
    end
  else
    % matching without wild card

    for iatom = 1:natom
      if strcmpi(strtrim(namelist(iatom, :)), query)
        logical_index(iatom) = true;
      end
    end
    
    % logical_index2 = strcmpi(strtrim(cellstr(namelist)), query);
    % logical_index = logical_index | logical_index2;

    % f = @trim_and_compare;
    % logical_index2 = cellfun(f, cellstr(namelist), cellstr(repmat(query, natom, 1)));
    % logical_index = logical_index | logical_index2;
  end
end


function s = trim_and_compare(name, query)
name = strtrim(name);
s = strcmpi(name, query);

