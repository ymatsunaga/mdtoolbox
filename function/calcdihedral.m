function s = calcdihedral(x);
%% calcdihedral
% calculate dihedral angles from Cartesian coordinates
%
%% Syntax
%# dihed = calcdihedral(x);
%
%% Description
% Calculate dihedral angles from the input Cartesian coordinates 
% of four atoms. 
% Input coordinates variable x has have 'nstep' rows and '3*4=12' columns.
% Each row has the XYZ coordinates of four atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3) x(4) y(4) z(4)].
%
% * x       - coordinates of four atoms [nstepx12 double]
% * dihed   - dihedral angles in degree [nstepx1 double]
%
%% Example
%# s = dihedral(x(:,1:(3*4)));
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
  s(istep) = acos(dot(m1,m2)./(norm(m1).*norm(m2))).*180./pi;
  rotdirection = dot(d2,cross(m1,m2));
  if rotdirection < 0
    s(istep) = - s(istep);
  end
end

