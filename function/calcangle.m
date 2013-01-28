function s = calcangle(x);
%% calcangle
% calculate the angles of three atoms from their Cartesian coordinates
%
%% Syntax
%# angle = calcangle(x);
%
%% Description
% Calculate angles from the input Cartesian coordinates 
% of three atoms in radians.
% Input coordinates variable x has have 'nstep' rows and '3*3=9' columns.
% Each row has the XYZ coordinates of three atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3)].
%
% * x       - coordinates of three atoms [nstep x 9 double]
% * angle   - angles in radians [nstep x 1 double]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# index = [5 7 9];
%# index3 = to3(index);
%# angle = calcangle(trj(:, index3)).*180./pi;
%
%% See alo
% calcbond, calcdihedral
% 

nstep = size(x,1);
s = zeros(nstep,1);

for istep = 1:nstep
  d1 = x(istep,1:3) - x(istep,4:6);
  d2 = x(istep,7:9) - x(istep,4:6);
  s(istep) = acos(dot(d1,d2)./(norm(d1).*norm(d2)));
end

