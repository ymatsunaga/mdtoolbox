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
% * trj        - coordinates of atoms [nframe x natom3]
% * quadruplet - quadruplet indices whose angles are calculated [nquadruplet x 3]
% * dihedral   - dihedral angles of the quadruplet (in radian) [nframe x nquadruplet]
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
if ~exist('quadruplet', 'var') || isempty(quadruplet)
  quadruplet = [1 2 3 4];
end

nframe = size(trj, 1);
nquadruplet = size(quadruplet, 1);

%% calculation
dihedral = zeros(nframe, nquadruplet);

for iquadruplet = 1:nquadruplet
  index1 = to3(quadruplet(iquadruplet, 1));
  index2 = to3(quadruplet(iquadruplet, 2));
  index3 = to3(quadruplet(iquadruplet, 3));
  index4 = to3(quadruplet(iquadruplet, 4));
  for iframe = 1:nframe
    d1 = trj(iframe, index1) - trj(iframe, index2);
    d2 = trj(iframe, index3) - trj(iframe, index2);
    d3 = trj(iframe, index3) - trj(iframe, index4);
    m1 = cross(d1, d2);
    m2 = cross(d2, d3);
    dihedral(iframe, iquadruplet) = acos(dot(m1, m2)./(norm(m1).*norm(m2)));
    rotdirection = dot(d2,cross(m1, m2));
    if rotdirection < 0
      dihedral(iframe, iquadruplet) = - dihedral(iframe, iquadruplet);
    end
  end
end

