function s = calcbond(x)
%% calcbond
% calculate the distances between two atom from their Cartesian coordinates
%
%% Syntax
%# bond = calcbond(x);
%
%% Description
% Calculate distances from the input Cartesian coordinates 
% of two atoms.
% Input coordinates variable x has have 'nstep' rows and '3*2=6' columns.
% Each row has the XYZ coordinates of three atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2)].
%
% * x      - coordinates of two atoms [nstep x 6 double]
% * bond   - distances between the two atoms [nstep x 1 double]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# index = [5 7];
%# index3 = to3(index);
%# bond = calcbond(trj(:, index3));
%
%% See alo
% calcangle, calcdihedral
% 

s = sqrt(sum((x(:, 1:3) - x(:, 4:6)).^2, 2));

