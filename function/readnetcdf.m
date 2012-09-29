function [trj, box, attributes] = readnetcdf(filename, index_atom, index_time)
%% readnetcdf
% read amber netcdf file
%
%% Syntax
%# trj = readnetcdf(filename);
%# trj = readnetcdf(filename, index_atom);
%# trj = readnetcdf(filename, [], index_time);
%# [trj, box] = readnetcdf(filename, index_atom, index_time);
%# [trj, box, attributes] = readnetcdf(filename, index_atom, index_time);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nstep' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input netcdf trajectory filename
% * index_atom - atom index or logical index specifying atoms to be read
% * index_time - time-step index or logical index specifying steps to be read
%                please note that this should be a regularly-spaced index
% * trj        - trajectory [nstep x natom3 double]
% * box        - box size [nstep x 3 double]
% * attributes - attributes (like header information) of netcdf file
%                [structure]
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
% support for velocities
% 

%% read attributes
finfo = ncinfo(filename);
attributes = struct;

for i = 1:numel(finfo.Attributes)
  attributes = setfield(attributes, finfo.Attributes(i).Name, finfo.Attributes(i).Value);
end
attributes

%% read dimensions
dimensions = struct;

for i = 1:numel(finfo.Dimensions)
  dimensions = setfield(dimensions, finfo.Dimensions(i).Name, finfo.Dimensions(i).Length);
end
dimensions

%% initialization
natom = dimensions.atom;
nstep = dimensions.frame;

if (nargin < 2) | (numel(index_atom) == 0)
  index_atom = 1:natom;
else
  if islogical(index_atom)
    index_atom = find(index_atom);
  end
  index_atom(index_atom > natom) = [];
end

if (nargin < 3) | (numel(index_time) == 0)
  index_time = 1:nstep;
else
  if islogical(index_time)
    index_time = find(index_time);
  end
  index_time(index_time > nstep) = [];
end

start_atom  = min(index_atom);
count_atom  = max(index_atom) - start_atom + 1;
stride_atom = 1;
index_atom = index_atom - start_atom + 1;

start_time  = min(index_time);
count_time  = numel(index_time);
stride_time = unique(diff(index_time));

%% read data
% trajectory in Angstrom
trj = ncread(filename, 'coordinates', [1 start_atom start_time], ...
             [3 count_atom count_time], [1 stride_atom stride_time]); 
trj = trj(:, index_atom, :);
[nspatial, natom, nstep] = size(trj);
trj = reshape(trj, nspatial*natom, nstep)';
trj = double(trj);

% box-size in Angstrom
if (nargout > 1)
  box = ncread(filename, 'cell_lengths', [1 start_time], ...
               [3 count_time], [1 stride_time]); 
  box = box';
end

% time in picosecond
% attribute.time = ncread(filename, 'time');

% box-angle in degree
% tmp = ncread(filename, 'cell_angles');
% attribute.angular = tmp';

