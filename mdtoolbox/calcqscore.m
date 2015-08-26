function qscore = calcqscore(trj, ref, id_residue, index_atom)
%% calcqscore
% calculate Q-score from given heavy atom coordinates
%
%% Syntax
%# qscore = calcqscore(trj, ref, id_residue, index_atom)
%
%% Description
%
% trj        -  input trajectory
% ref        -  coordinate of native structure
% id_residue - residue ID for atoms
% index_atom - index for heavy atoms. If omitted, all atom atoms in trj or ref are assumed to be heavy atoms
% 
% qscore     - Q-score [nframe x 1 double]
%
%% Example
%# [pdb, ref] = readpdb('run.pdb');
%# index_atom = selectid(pdb.resseq, 10:30) & ~selectname(pdb.name, 'H*');
%# trj = readnetcdf('run.nc');
%# q = calcqscore(trj, ref, pdb.resseq, index_atom);
% 
%% See also
%
%% References
% Definition of Q-score from heavy atoms is based on 
% [1] R. B. Best, G. Hummer, and W. A. Eaton, PNAS 110, 17874 (2013).
%

%% setup
nframe = size(trj, 1);

if exist('index_atom', 'var')
  trj = trj(:, to3(index_atom));
  ref = ref(to3(index_atom));
  id_residue = id_residue(index_atom);
end  

natom3 = size(trj, 2);
natom = natom3/3;
assert(natom3 == numel(ref), 'sizes of trj and ref are not consistent.');
assert(natom == numel(id_residue), 'sizes of trj and id_residue are not consistent.');

%% find native contacts for all heavy atoms
pair = [];
dist0 = [];

pair_all = nchoosek(1:natom, 2);
ires = id_residue(pair_all(:, 1));
jres = id_residue(pair_all(:, 2));
id = abs(ires - jres) > 3;
pair = pair_all(id, :);

dist0 = calcbond(ref, pair);
id = (dist0 < 4.5);
dist0 = dist0(id);
pair = pair(id, :);
disp(sprintf('Message: %d native contacts found', sum(id)));

%% calc Q-score
dist = calcbond(trj, pair);

BETA = 5.0;
LAMBDA = 1.8;
qscore = mean(1.0./(1.0 + exp(BETA * bsxfun(@minus, dist, LAMBDA*dist0))), 2);

