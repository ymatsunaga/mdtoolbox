function logical_index = selectid(id, index)
%% selectid
% Used for atom selection. Returns logical indices of atoms or residues which matched given IDs.
%
%% Syntax
%# logical_index = selectid(id, index)
%
%% Description
% Used for atom selection. 
% Compares atom or residue 'id's with index and returns
% logical index which have 'trues' for the matched atoms.
%
% * id            - atom or residue ids [natom x 1 integer]
% * index         - index [n x 1 (or 1 x n) integer]
% * logical_index - logical values which has trues for matched atoms
%                   [1 x natom logical]
%
%% Example
%# % read amber parm file
%# parm = readparm('ak.parm');
%#
%# % choose the atom of 1-3 residues
%# logical_index = selectid(parm.residue_id, 1:3);
%# index = find(logical_index);
%

% calc matches between id and index
logical_index = false(size(id));

for i = 1:numel(index)
  logical_index = logical_index | (id == index(i));
end

