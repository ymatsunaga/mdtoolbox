function contactMap = calccontactmap(crd, cutoff)
%% calccontactmap
%
%% Syntax
%#
%
%% Description
%
%% Example
%# example: 
%# crd = readpdb('ak.pdb');
%# contactMap = calccontactmap(crd, 8.0)
%# imagesc(contactMap);
%# axis xy;
%# colorbar
% 
%% See also
%
%% References
%

distanceMatrix = calcdistancematrix(crd);
contactMap = double(distanceMatrix < cutoff);

