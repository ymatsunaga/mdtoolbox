function [trj, box, vel, temp, attributes] = readnetcdf(filename, index_atom, index_time)
%% readnetcdf
% read amber netcdf file
%
%% Syntax
%# trj = readnetcdf(filename);
%# trj = readnetcdf(filename, index_atom);
%# trj = readnetcdf(filename, [], index_time);
%# [trj, box] = readnetcdf(filename, index_atom, index_time);
%# [trj, box, vel] = readnetcdf(filename, index_atom, index_time);
%# [trj, box, vel, temp] = readnetcdf(filename, index_atom, index_time);
%# [trj, box, vel, attributes] = readnetcdf(filename, index_atom, index_time);
%
%% Description
% The XYZ coordinates of atoms are read into 'trj' variable
% which has 'nframe' rows and '3*natom' columns.
% Each row of 'trj' has the XYZ coordinates of atoms in order 
% [x(1) y(1) z(1) x(2) y(2) z(2) ... x(natom) y(natom) z(natom)].
%
% * filename   - input netcdf trajectory filename
% * index_atom - atom index or logical index specifying atoms to be read
% * index_time - time-frame index or logical index specifying frames to be read
%                please note that this should be a regularly-spaced index
% * trj        - coordinates in units of angstrom [nframe x natom3 double]
% * box        - box size in units of angstrom [nframe x 3 double]
% * vel        - velocity in units of angstrom/picosecond [nframe x natom3 double]
% * temp       - target temperature in units of kelvin [nframe x 1 double]
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
%   double cell_lengths(frame, cell_spatial) units="angstrom"
%   double cell_angles(frame, cell_angular) units="degree"
%   float velocities(frame, atom, spatial) units="angstrom/picosecond"
%                                           scale_factor=20.455f
%% TODO
%
% 

%% displays groups, dimensions, variable definitions, and all attributes in the NetCDF data
%ncdisp(filename);

%% read attributes
finfo = ncinfo(filename);
attributes = struct;

for i = 1:numel(finfo.Attributes)
  attributes = setfield(attributes, finfo.Attributes(i).Name, finfo.Attributes(i).Value);
end

%% read dimensions
dimensions = struct;

for i = 1:numel(finfo.Dimensions)
  dimensions = setfield(dimensions, finfo.Dimensions(i).Name, finfo.Dimensions(i).Length);
end

%% determine flags
is_trj = false;
is_box = false;
is_vel = false;
is_temp = false;

for i = 1:numel(finfo.Variables)
  varname = finfo.Variables(i).Name;
  if strncmpi(varname, 'coordinates', numel('coordinates'))
    is_trj = true;
  end
  if strncmpi(varname, 'cell_lengths', numel('cell_lengths'))
    is_box = true;
  end
  if strncmpi(varname, 'velocities', numel('velocities'))
    is_vel = true;
  end
  if strncmpi(varname, 'temp0', numel('temp0'))
    is_temp = true;
  end
end

if nargout < 2
  is_box = false;
end

if nargout < 3
  is_vel = false;
end

if nargout < 4
  is_temp = false;
end

%% initialization
natom = dimensions.atom;
nframe = dimensions.frame;

if ~exist('index_atom', 'var') || isempty(index_atom)
  index_atom = 1:natom;
else
  if islogical(index_atom)
    index_atom = find(index_atom);
  end
  index_atom(index_atom > natom) = [];
end

if ~exist('index_time', 'var') || isempty(index_time)
  index_time = 1:nframe;
else
  if islogical(index_time)
    index_time = find(index_time);
  end
  index_time(index_time > nframe) = [];
end

start_atom  = min(index_atom);
count_atom  = max(index_atom) - start_atom + 1;
stride_atom = 1;
index_atom = index_atom - start_atom + 1;

start_time  = min(index_time);
count_time  = numel(index_time);
stride_time = unique(diff(index_time));

%% read data
% coordinates in Angstrom
if is_trj
  trj = ncread(filename, 'coordinates', [1 start_atom start_time], [3 count_atom count_time], [1 stride_atom stride_time]); 
  trj = trj(:, index_atom, :);
  [nspatial, natom, nframe] = size(trj);
  trj = reshape(trj, nspatial*natom, nframe)';
  trj = double(trj);
else
  trj = [];
end

% box-sizes in Angstrom
if is_box
  box = ncread(filename, 'cell_lengths', [1 start_time], [3 count_time], [1 stride_time]); 
  box = box';
else
  box = [];
end

% velocities in Angstrom/picosecond
if is_vel
  vel = ncread(filename, 'velocities', [1 start_atom start_time], [3 count_atom count_time], [1 stride_atom stride_time]); 
  vel = vel(:, index_atom, :);
  [nspatial, natom, nframe] = size(vel);
  vel = reshape(vel, nspatial*natom, nframe)';
  vel = double(vel);
  vel = 20.455*vel;
else
  vel = [];
end

% target temperatures in Kelvin
if is_temp
  temp = ncread(filename, 'temp0', [start_time], [count_time], [stride_time]); 
else
  temp = [];
end

