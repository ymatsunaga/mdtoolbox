function [pair, dist] = calccontact(crd, atom_name, residue_id, cutoff, sub_residue_id)
%% calccontact
% calculate contact residue pairs from an all-atom structure
%
%% Syntax
%# [pair, dist] = calccontact(crd, atom_name, residue_id, cutoff)
%# [pair, dist] = calccontact(crd, atom_name, residue_id, cutoff, sub_residue_id)
%
%% Description
% A residue pair is defined as “contacted” when at least one
% non-hydrogen atom in residue i is within 'cutoff' A from a
% non-hydrogen atom of residue j.
% If 'sub_residue_id' is given, contacts are searched only for A
% and B residues which belong to 'sub_residue_id'.
%
%% Example
%# [pdb, crd] = readpdb('lys.pdb');
%# index = selectname(pdb.record, 'ATOM');
%# pdb = substruct(pdb, index);
%# crd = crd(to3(index));
%# [pair, dist] = calccontact(crd, pdb.name, pdb.resseq);
%# spy(sparse(pair(:,1), pair(:,2), 1, 164, 164))
% 
%% See also
%
%% References
%

%% setup
if nargin < 4
  cutoff = 6.5;
end

if nargin > 4
  residue_unique = unique(sort(sub_residue_id));
else
  residue_unique = unique(sort(residue_id));
end

% choose non-hydrogen atoms
index = ~selectname(atom_name, 'H*');
index = find(index);

crd          = crd(to3(index));
atom_name    = atom_name(index, :);
residue_id   = residue_id(index);

%% search contact residue pairs
pair = [];
dist = [];
for i = 1:numel(residue_unique)
  ires = residue_unique(i);
  index = selectrange(crd, residue_id == ires, cutoff);
  residue_contact = unique(residue_id(index));
  for j = 1:numel(residue_contact)
    jres = residue_contact(j);
    if jres > (ires + 3)
      pair = [pair; [ires jres]];
      %Ca-Ca distance
      if nargout >= 2
        ca_i = find(selectid(residue_id, ires) & selectname(atom_name, 'CA'));
        ca_j = find(selectid(residue_id, jres) & selectname(atom_name, 'CA'));
        assert(numel(ca_i) == 1, sprintf('failed to identify CA of %d-th residue', ires));
        assert(numel(ca_j) == 1, sprintf('failed to identify CA of %d-th residue', jres));
        d = sqrt(sum((crd(to3(ca_i)) - crd(to3(ca_j))).^2));
        dist = [dist; d];
      end
    end
  end
end

