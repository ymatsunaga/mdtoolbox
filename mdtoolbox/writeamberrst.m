function writeamberrst(filename, crd, box, vel, title)
%% writenetcdf
% write amber netcdf restart file
%
%% Syntax
%# writeamberrst(filename, crd, box);
%# writeamberrst(filename, crd, box, vel);
%# writeamberrst(filename, crd, box, vel, title);
%
%% Description
% This code puts coordinates and velocities into an amber netcdf restart file. 
%
% * filename  - output filename
% * crd       - coordinates [1 x natom3 double or single]
% * box       - box size [1 x 3 double]
% * crd       - velocities [1 x natom3 double or single]
% * title     - title [chars]
%
%% Example
% crd(:, 1:3:end) = crd(:, 1:3:end) + 10;
% writeamberrst('ak_translated.nc', crd, box, vel);
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
%   Conventions (required): should be "AMBERRESTART"
%   ConventionVersion (required): "1.0" or later
%   application (optional): ex. "AMBER"
%   program (optional): ex. "sander"
%   programVersion (optional): ex. "10.0"
%   title (optional): ex. "productino run"
% Dimensions:
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
%   double time units="picosecond"
%   double coordinates(atom, spatial) units="angstrom"
%   double cell_length(cell_spatial) units="angstrom"
%   double cell_angles(cell_angular) units="degree"
%   double velocities(atom, spatial) units="angstrom/picosecond"
%                                    scale_factor=20.455f
% 
%% TODO
% remd support, etc.
%

%% initialization
scale_factor = 20.455;
natom3 = numel(crd);
natom = natom3 / 3;

if ~iscolumn(crd)
  crd = crd';
end
crd = reshape(crd, 3, natom);

if ~iscolumn(vel)
  vel = vel';
end
vel = reshape(vel, 3, natom);

if ~exist('title', 'var')
  title = 'CREATED BY MATLAB';
end

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% open file
cmode1 = netcdf.getConstant('CLOBBER');
cmode2 = netcdf.getConstant('64BIT_OFFSET');
cmode = bitor(cmode1, cmode2);
ncid = netcdf.create(filename, cmode);

%% global attributes
varid_global = netcdf.getConstant('GLOBAL');

netcdf.putAtt(ncid, varid_global, 'ConventionVersion', '1.0');
netcdf.putAtt(ncid, varid_global, 'Conventions', 'AMBERRESTART');
netcdf.putAtt(ncid, varid_global, 'application', 'MDTOOLBOX');
netcdf.putAtt(ncid, varid_global, 'program', 'MATLAB');
netcdf.putAtt(ncid, varid_global, 'programVersion', '1.0');

if ~exist('title', 'var')
  title = 'CREATED BY MATLAB';
end
netcdf.putAtt(ncid, varid_global, 'title', title);

%% dimensions
%dimid_frame          = netcdf.defDim(ncid, 'frame', netcdf.getConstant('NC_UNLIMITED'));
%dimid_remd_dimension = netcdf.defDim(ncid, 'remd_dimension', 2);
dimid_spatial         = netcdf.defDim(ncid, 'spatial', 3);
dimid_atom            = netcdf.defDim(ncid, 'atom', natom);
dimid_cell_spatial    = netcdf.defDim(ncid, 'cell_spatial', 3);
dimid_cell_angular    = netcdf.defDim(ncid, 'cell_angular', 3);
dimid_label           = netcdf.defDim(ncid, 'label', 5);

%% variables
varid_spatial      = netcdf.defVar(ncid, 'spatial', netcdf.getConstant('NC_CHAR'), dimid_spatial);
varid_cell_spatial = netcdf.defVar(ncid, 'cell_spatial', netcdf.getConstant('NC_CHAR'), dimid_cell_spatial);
varid_cell_angular = netcdf.defVar(ncid, 'cell_angular', netcdf.getConstant('NC_CHAR'), [dimid_label, dimid_cell_angular]);

varid_time         = netcdf.defVar(ncid, 'time', netcdf.getConstant('NC_DOUBLE'), []);
netcdf.putAtt(ncid, varid_time, 'units', 'picosecond');

varid_coordinates  = netcdf.defVar(ncid, 'coordinates', netcdf.getConstant('NC_DOUBLE'), [dimid_spatial, dimid_atom]);
netcdf.putAtt(ncid, varid_coordinates, 'units', 'angstrom');

varid_cell_lengths = netcdf.defVar(ncid, 'cell_lengths', netcdf.getConstant('NC_DOUBLE'), dimid_cell_spatial);
netcdf.putAtt(ncid, varid_cell_lengths, 'units', 'angstrom');

varid_cell_angles  = netcdf.defVar(ncid, 'cell_angles', netcdf.getConstant('NC_DOUBLE'), dimid_cell_angular);
netcdf.putAtt(ncid, varid_cell_angles, 'units', 'degree');

varid_velocities   = netcdf.defVar(ncid, 'velocities', netcdf.getConstant('NC_DOUBLE'), [dimid_spatial, dimid_atom]);
netcdf.putAtt(ncid, varid_velocities, 'units', 'angstrom/picosecond');
netcdf.putAtt(ncid, varid_velocities, 'scale_factor', scale_factor);

%varid_forces       = netcdf.defVar(ncid, 'forces', netcdf.getConstant('NC_DOUBLE'), [dimid_atom, dimid_spatial]);
%netcdf.putAtt(ncid, varid_forces, 'units', 'amu*angstrom/picosecond^2');

%varid_temp0        = netcdf.defVar(ncid, 'temp0', netcdf.getConstant('NC_DOUBLE'), []);
%netcdf.putAtt(ncid, varid_tmp0, 'units', 'kelvin');

%varid_remd_type    = netcdf.defVar(ncid, 'remd_type', netcdf.getConstant('NC_INT'), []);
%varid_remd_indices = netcdf.defVar(ncid, 'remd_indices', netcdf.getConstant('NC_INT'), []);

netcdf.endDef(ncid);

%% wtite data
netcdf.putVar(ncid, varid_spatial, ['xyz']);
netcdf.putVar(ncid, varid_time, 0.0); % TODO
netcdf.putVar(ncid, varid_coordinates, crd);

if exist('vel', 'var') && ~isempty(vel)
  netcdf.putVar(ncid, varid_velocities, vel./scale_factor);
end

if exist('box', 'var') && ~isempty(box)
  netcdf.putVar(ncid, varid_cell_spatial, ['abc']);
  netcdf.putVar(ncid, varid_cell_angular, ['alpha'; 'beta '; 'gamma']');
  netcdf.putVar(ncid, varid_cell_lengths, box');
  netcdf.putVar(ncid, varid_cell_angles, [90.0 90.0 90.0]');
end

%% close file
netcdf.close(ncid);

