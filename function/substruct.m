function str_sub = substruct(str, index)
%% substruct
% create a subset structure from a structure of arrays of same size
%
%% Syntax
%# str_sub = substruct(str, index)
%
%% Description
% This code assumes that the given structure has
% arrays of same size, and just create a subset 
% of these arrays
%
% * str     - structure
% * index   - index or logical index
% * str_sub - structure
%
%% Example
%# pdb = readpdb('ak.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb_ca = substruct(pdb, index_ca);
% 

str_sub = structfun(@(x) x(index, :), str, 'uniformOutput', false);

