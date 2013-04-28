function s = calcdihedral(x);
%% calcdihedral
% calculate the dihedral angles of four atoms from their Cartesian coordinates
%
%% Syntax
%# dihed = calcdihedral(x);
%
%% Description
% Calculate dihedral angles from the input Cartesian coordinates 
% of four atoms in radians.
% Input coordinates variable x has have 'nstep' rows and '3*4=12' columns.
% Each row has the XYZ coordinates of four atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3) x(4) y(4) z(4)].
%
% * x       - coordinates of four atoms [nstep x 12 double]
% * dihed   - dihedral angles in radians [nstep x 1 double]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# index_phi = [5 7 9 15];
%# index_psi = [7 9 15 17];
%# index_phi3 = to3(index_phi);
%# index_psi3 = to3(index_psi);
%# phi = calcdihedral(trj(:, index_phi3)).*180./pi;
%# psi = calcdihedral(trj(:, index_psi3)).*180./pi;
%
%% See alo
% calcbond, calcangle
% 

nstep = size(x,1);
s = zeros(nstep,1);

for istep = 1:nstep
  d1 = x(istep,1:3) - x(istep,4:6);
  d2 = x(istep,7:9) - x(istep,4:6);
  d3 = x(istep,7:9) - x(istep,10:12);
  m1 = cross(d1,d2);
  m2 = cross(d2,d3);
  s(istep) = acos(dot(m1,m2)./(norm(m1).*norm(m2)));
  rotdirection = dot(d2,cross(m1,m2));
  if rotdirection < 0
    s(istep) = - s(istep);
  end
end

