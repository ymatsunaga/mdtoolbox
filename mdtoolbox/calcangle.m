function angle = calcangle(trj, triplet)
%% calcangle
% calculate angles of triplet atoms from their Cartesian coordinates
%
%% Syntax
%# angle = calcangle(trj);
%# angle = calcangle(trj, triplet);
%
%% Description
% Calculate angles from the input trajectory of Cartesian coordinates.
% Triplets, whose angles are calculated, can be specified via the
% variable (triplet).
%
% * trj     - coordinates of atoms [nframe x natom3]
% * triplet - triplet indices whose angles are calculated [ntriplet x 3]
% * angle   - angles of the triplets (in radian) [nframe x ntriplet]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# triplet = [5 7 9; 7 9 10];
%# angle = calcangle(trj, triplet).*180./pi;
%
%% See alo
% calcbond, calcdihedral
% 

%% initialization
if ~exist('triplet', 'var') || isempty(triplet)
  triplet = [1 2 3];
end

nframe = size(trj, 1);
ntriplet = size(triplet, 1);

%% calculation
angle = zeros(nframe, ntriplet);

for itriplet = 1:ntriplet
  index1 = to3(triplet(itriplet, 1));
  index2 = to3(triplet(itriplet, 2));
  index3 = to3(triplet(itriplet, 3));
  for iframe = 1:nframe
    d1 = trj(iframe, index1) - trj(iframe, index2);
    d2 = trj(iframe, index3) - trj(iframe, index2);
    angle(iframe, itriplet) = acos(dot(d1, d2)./(norm(d1).*norm(d2)));
  end
end

