function str_union = addstruct(str1, str2)
%% addstruct
% create a structure by making the union of arrays of two structure variables
%
%% Syntax
%# str_union = addstruct(str1, str2)
%
%% Description
% This routine assumes that the given individual structure (str1, or
% str2) has arrays of same size. The sizes can be different in str1
% and str2. 
%
% * str1      - structure
% * str2      - structure
% * str_union - structure
%
%% Example
%# pdb1 = readpdb('enzyme.pdb');
%# pdb2 = readpdb('substrate.pdb');
%# pdb_complex = addstruct(pdb1, pdb2);
% 

str_union = struct;
s = fieldnames(str1);
for i = 1:numel(s)
  c = s{i};
  str_union.(c) = [str1.(c); str2.(c)];
end

