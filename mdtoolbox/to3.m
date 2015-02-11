function index3 = to3(index)
%% to3
% convert 1...N atom-index to 1...3N xyz-index
%
%% Syntax
%# index3 = to3(index)
%
%% Description
% This code converts 1...N atom-index to 1...3N xyz-index.
% Useful for the operation of Cartesian coordinates. 
%
% * index  - array of integers 
%            (Ex. [1 3])
%            or
%            logcal indexing
%            (Ex. logical([1 0 1]))
% * index3 - array of integers 
%            (Ex. [1 2 3 7 8 9])
%            or
%            logical indexing
%            (Ex. logical[1 1 1 0 0 0 1 1 1])
%
%% Example
%# pdb = readpdb('ak.pdb');
%# index = selectname(pdb.name, 'CA');
%# index3 = to3(index);
%# trj = trj(:,index3);
% 

if islogical(index)
  index3 = false(1, numel(index)*3);
  index3(1:3:end) = index;
  index3(2:3:end) = index;
  index3(3:3:end) = index;

else
  index3 = zeros(1, numel(index)*3);
  index3(1:3:end) = 3.*(index-1) + 1;
  index3(2:3:end) = 3.*(index-1) + 2;
  index3(3:3:end) = 3.*(index-1) + 3;

end

