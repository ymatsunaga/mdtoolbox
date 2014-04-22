function dihedral = calcdihedral(trj, quadruplet);
%% calcdihedral
% calculate dihedral angles of quadruplet atoms from their Cartesian coordinates
%
%% Syntax
%# dihedral = calcdihedral(trj);
%# dihedral = calcdihedral(trj, quadruplet);
%
%% Description
% Calculate dihedral angles from the input trajectory of Cartesian coordinates.
% Quadruplets, whose dihedral angles are calculated, can be specified via the
% variable (quadruplet).
%
% * trj        - coordinates of atoms [nstep x natom3]
% * quadruplet - quadruplet indices whose angles are calculated [nquadruplet x 3]
% * dihedral   - dihedral angles of the quadruplet (in radian) [nstep x nquadruplet]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# quadruplet = [5 6 7 8; 6 7 8 9];
%# dihedral = calcdihedral(trj, quadruplet).*180./pi;
%
%% See alo
% calcbond, calcangle
% 

%% initialization
if ~exist('quadruplet', 'var')
  triplet = [1 2 3 4];
end

nstep = size(trj, 1);
nquadruplet = size(quadruplet, 1);

%% calculation
dihedral = zeros(nstep, nquadruplet);

for iquadruplet = 1:nquadruplet
  index1 = to3(quadruplet(iquadruplet, 1));
  index2 = to3(quadruplet(iquadruplet, 2));
  index3 = to3(quadruplet(iquadruplet, 3));
  index4 = to3(quadruplet(iquadruplet, 4));
  for istep = 1:nstep
    d1 = trj(istep, index1) - trj(istep, index2);
    d2 = trj(istep, index3) - trj(istep, index2);
    d3 = trj(istep, index3) - trj(istep, index4);
    m1 = cross(d1, d2);
    m2 = cross(d2, d3);
    dihedral(istep, iquadruplet) = acos(dot(m1, m2)./(norm(m1).*norm(m2)));
    rotdirection = dot(d2,cross(m1, m2));
    if rotdirection < 0
      dihedral(istep, iquadruplet) = - dihedral(istep, iquadruplet);
    end
  end
end

