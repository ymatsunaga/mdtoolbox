function rc = writenetcdf(filename, trj, box, title)
%% writenetcdf
% write amber netcdf file
%
%% Syntax
%# writenetcdf(filename, trj);
%# writenetcdf(filename, trj, box);
%# writenetcdf(filename, trj, box, title);
%# writenetcdf(filename, trj, [], title);
%
%% Description
% This code puts trajectories into a dcd file. 
% If header information is not given, 
% default values are assumed. 
%
% * filename  - output dcd trajectory filename
% * trj       - trajectory [nstep x natom3 double]
% * box       - box size [nstep x 3 double]
% * title     - title [chars]
%
%% Example
% trj = readnetcdf('ak.nc');
% trj(:, 1:3:end) = trj(:, 1:3:end) + 10;
% writenetcdf('ak_translated.nc', trj);
% 
%% References
% For the file format specifications, see
% http://ambermd.org/netcdf/nctraj.html
%
% Here is a brief summary of the format:
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
%   char spatial(spatial): "xyz", usually
%   char cell_spatial(cell_spatial): "abc", usually
%   char cell_angular(cell_angular, label): "alpha", "beta ", "gamma", usually
%   data variables:
%   float time(frame) units="picosecond"
%   float coordinates(frame, atom, spatial) units="angstrom"
%   double cell_length(frame, cell_spatial) units="angstrom"
%   double cell_angles(frame, cell_angular) units="degree"
%   double velocities(frame, atom, spatial) units="angstrom/picosecond"
%                                           scale_factor=20.455f
% 

[nStep, nAtom3] = size(trj);
nAtom = nAtom3 / 3;

trj = single(trj);
trj = trj';
trj = reshape(trj,3,nAtom,nStep);

% create 64 bit offset netcdf file
ncID = netcdf.create(filename, '64BIT_OFFSET');

% open the file
%ncID = netcdf.open(filename, 'NC_WRITE');

% define dimensions
%dimID_frame        = netcdf.defDim(ncID, 'frame', netcdf.getConstant('NC_UNLIMITED'));
dimID_frame        = netcdf.defDim(ncID, 'frame', nStep);
dimID_spatial      = netcdf.defDim(ncID, 'spatial', 3);
dimID_atom         = netcdf.defDim(ncID, 'atom', nAtom);
if nargin >= 3
  dimID_cell_spatial = netcdf.defDim(ncID, 'cell_spatial', 3);
  dimID_cell_angular = netcdf.defDim(ncID, 'cell_angular', 3);
  dimID_label        = netcdf.defDim(ncID, 'label', 5);
end

% define new variables
%varID              = netcdf.defVar(ncID,'my_var','double', dimID);
varID_spatial      = netcdf.defVar(ncID, 'spatial', 'char', dimID_spatial);
varID_time         = netcdf.defVar(ncID, 'time', 'float', dimID_frame);
varID_coordinates  = netcdf.defVar(ncID, 'coordinates', 'float', [dimID_spatial dimID_atom dimID_frame]);
if nargin >= 3
  varID_cell_spatial = netcdf.defVar(ncID, 'cell_spatial', 'char', dimID_cell_spatial);
  varID_cell_angular = netcdf.defVar(ncID, 'cell_angular', 'char', [dimID_label dimID_cell_angular]);
  varID_cell_lengths = netcdf.defVar(ncID, 'cell_lengths', 'double', [dimID_spatial dimID_frame]);
  varID_cell_angles  = netcdf.defVar(ncID, 'cell_angles', 'double', [dimID_spatial dimID_frame]);
end

% leave define mode and enter data mode to write data
netcdf.endDef(ncID);

% write data to variable.
netcdf.putVar(ncID, varID_spatial, ['xyz']');
netcdf.putVar(ncID, varID_time, 0, nStep, [1:nStep]);
netcdf.putVar(ncID, varID_coordinates, [0 0 0], [3 nAtom nStep], trj);
if nargin >= 3
  netcdf.putVar(ncID, varID_cell_spatial, ['abc']');
  netcdf.putVar(ncID, varID_cell_angular, ['alpha'; 'beta '; 'gamma']');
  netcdf.putVar(ncID, varID_cell_lengths, [0 0], [3 nStep], box');
  netcdf.putVar(ncID, varID_cell_angles, [0 0], [3 nStep], repmat(90,3,nStep));
end

% Re-enter define mode.
netcdf.reDef(ncID);

% Create attributes associated with the variables.
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'title', title);
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'application', 'MATLAB');
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'program', 'writenetcdf()');
%netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'programVersion', '10.0');
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'programVersion', version);
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'AMBER');
netcdf.putAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'ConventionVersion', '1.0');

netcdf.close(ncID);


