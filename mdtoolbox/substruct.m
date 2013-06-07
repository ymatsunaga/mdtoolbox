function str_sub = substruct(str, index)
%% substruct
% create a subset structure from a structure of arrays
%
%% Syntax
%# str_sub = substruct(str, index)
%
%% Description
% This code assumes that the given structure variable has
% arrays of same size, and create a subset of these arrays
%
% * str     - structure
% * index   - index or logical index
% * str_sub - structure
%
%% Example
%# pdb = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb_ca = substruct(pdb, index_ca);
%# writepdb('lys_ca.pdb', pdb_ca);
% 

str_sub = structfun(@(x) x(index, :), str, 'UniformOutput', false);

