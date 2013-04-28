function distanceMatrix = calcdistancematrix(crd)
%% calcdistancematrix
%
%% Syntax
%#
%
%% Description
%
%% Example
%# crd = readpdb('ak.pdb');
%# distanceMatrix = calcdistancematrix(crd)
%# imagesc(distanceMatrix);
%# axis xy;
%# colorbar
% 
%% See also
%
%% References
% 

c1 = bsxfun(@minus,crd(1:3:end),crd(1:3:end)');
c2 = bsxfun(@minus,crd(2:3:end),crd(2:3:end)');
c3 = bsxfun(@minus,crd(3:3:end),crd(3:3:end)');

distanceMatrix = sqrt(c1.^2 + c2.^2 + c3.^2);

