function [trj, box, attribute] = readnetcdf(filename, index)
%% readnetcdf
% read amber netcdf file
%
%% Syntax
%# trj = readnetcdf(filename);
%# trj = readnetcdf(filename, index);
%# [trj, box] = readnetcdf(filename, index);
%# [trj, box, attribute] = readnetcdf(filename, index);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nstep' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename  - input netcdf trajectory filename
% * index     - index or logical index specifying atoms to be read
% * trj       - trajectory [nstep x natom3 double]
% * box       - box size [nstep x 3 double]
% * attribute - attributes (like header information) of netcdf file
%               [structure]
%
%% Example
%# trj = readnetcdf('ak.nc');
%
%% See also
% writenetcdf
%
%% References
% For the file format specifications, see
% http://ambermd.org/netcdf/nctraj.html
%
% Here is a summary of the format:
% Encoding: 
%   NetCDF version 3.x, not the HDF5 encoding.
%   64bit offsets
% Global attributes:
%   Conventions (required): should be "AMBER"
%   ConventionVersion (required): "1.0" or later
%   application (optional): ex. "AMBER"
%   program (optional): ex. "sander"
%   programVersion (optional): ex. "10.0"
%   title (optional): ex. "productino run"
% Dimensions:
%   frame (required, length unlimited): "UNLIMITED", usually
%   spatial (required, length 3): 3, usually
%   atom (required, length set as appropriate): # of atoms
%   cell_spatial (optional, length 3): 3, usually
%   cell_angular (optional, length 3): 3, usually
%   label (optional, length set as appropriate): length of the longest label string, usually 5
% Variables:
%   label variables:
%   char spatial(spatial): "xyz"
%   char cell_spatial(cell_spatial): "abc"
%   char cell_angular(cell_angular, label): "alpha", "beta ", "gamma"
%   data variables:
%   float time(frame) units="picosecond"
%   float coordinates(frame, atom, spatial) units="angstrom"
%   double cell_length(frame, cell_spatial) units="angstrom"
%   double cell_angles(frame, cell_angular) units="degree"
%   double velocities(frame, atom, spatial) units="angstrom/picosecond"
%                                           scale_factor=20.455f
%% TODO
% support for reading subset snapshots
% 

% read attributes
attribute = struct([]);
attributeNames = {'Conventions', ...
                   'ConventionVersion', ...
                   'application', ...
                   'program', ...
                   'programVersion', ...
                   'title'};
for name = attributeNames
  name = name{1};
  attribute(1).(name) = ncreadatt(filename, '/', name);
end
attribute

% trajectory in Angstrom
trj = ncread(filename, 'coordinates'); 
[nSpatial, nAtom, nStep] = size(trj);
trj = reshape(trj, nSpatial*nAtom, nStep)';
trj = double(trj);

% box-size in Angstrom
box = ncread(filename, 'cell_lengths'); 
box = box';

% time in picosecond
attribute.time = ncread(filename, 'time');

% box-angle in degree
tmp = ncread(filename, 'cell_angles');
attribute.angular = tmp';


