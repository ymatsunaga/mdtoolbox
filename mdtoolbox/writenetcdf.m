function writenetcdf(filename, trj, box, title)
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
% This code puts trajectories into a netcdf file. 
%
% * filename  - output dcd trajectory filename
% * trj       - trajectory [nframe x natom3 double or single]
% * box       - box size [nframe x 3 double]
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
%% TODO
% support for velocities
%

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% initialization
[nframe, natom3] = size(trj);
natom = natom3 / 3;

trj = trj';
trj = reshape(trj, 3, natom, nframe);

if exist('box', 'var') && ~isempty(box)
  if (size(box, 1) == 1)
    box = repmat(box, nframe, 1);
  end
end

if ~exist('title', 'var')
  title = 'CREATED BY MATLAB';
end

%% write netcdf schema
finfo.Name   = '/';
finfo.Format = '64bit';
finfo.Groups = [];

finfo.Attributes(1).Name  = 'title';
finfo.Attributes(1).Value = title;

finfo.Attributes(2).Name  = 'application';
finfo.Attributes(2).Value = 'MATLAB';

finfo.Attributes(3).Name  = 'program';
finfo.Attributes(3).Value = 'MATLAB';

finfo.Attributes(4).Name  = 'programVersion';
finfo.Attributes(4).Value = version;

finfo.Attributes(5).Name  = 'Conventions'; %required
finfo.Attributes(5).Value = 'AMBER';

finfo.Attributes(6).Name  = 'ConventionVersion';
finfo.Attributes(6).Value = '1.0';

finfo.Dimensions(1).Name      = 'frame';
finfo.Dimensions(1).Length    = nframe;
finfo.Dimensions(1).Unlimited = true;

finfo.Dimensions(2).Name      = 'spatial';
finfo.Dimensions(2).Length    = 3;
finfo.Dimensions(2).Unlimited = false;

finfo.Dimensions(3).Name      = 'atom';
finfo.Dimensions(3).Length    = natom;
finfo.Dimensions(3).Unlimited = false;

if exist('box', 'var') && ~isempty(box)
  finfo.Dimensions(4).Name      = 'label';
  finfo.Dimensions(4).Length    = 5;
  finfo.Dimensions(4).Unlimited = false;

  finfo.Dimensions(5).Name      = 'cell_spatial';
  finfo.Dimensions(5).Length    = 3;
  finfo.Dimensions(5).Unlimited = false;

  finfo.Dimensions(6).Name      = 'cell_angular';
  finfo.Dimensions(6).Length    = 3;
  finfo.Dimensions(6).Unlimited = false;
end

finfo.Variables(1).Name                 = 'spatial';
finfo.Variables(1).Dimensions.Name      = 'spatial';
finfo.Variables(1).Dimensions.Length    = 3;
finfo.Variables(1).Dimensions.Unlimited = false;
finfo.Variables(1).Size                 = 3;
finfo.Variables(1).Datatype             = 'char';
finfo.Variables(1).Attributes           = [];
finfo.Variables(1).ChunkSize            = [];
finfo.Variables(1).FillValue            = [];
finfo.Variables(1).DeflateLevel         = [];
finfo.Variables(1).Shuffle              = false;

finfo.Variables(2).Name                 = 'time';
finfo.Variables(2).Dimensions.Name      = 'frame';
finfo.Variables(2).Dimensions.Length    = nframe;
finfo.Variables(2).Dimensions.Unlimited = true;
finfo.Variables(2).Size                 = nframe;
finfo.Variables(2).Datatype             = 'single';
finfo.Variables(2).Attributes.Name      = 'units';
finfo.Variables(2).Attributes.Value     = 'picosecond';
finfo.Variables(2).ChunkSize            = [];
finfo.Variables(2).FillValue            = [];
finfo.Variables(2).DeflateLevel         = [];
finfo.Variables(2).Shuffle              = false;

finfo.Variables(3).Name                    = 'coordinates';
finfo.Variables(3).Dimensions(1).Name      = 'spatial';
finfo.Variables(3).Dimensions(1).Length    = 3;
finfo.Variables(3).Dimensions(1).Unlimited = false;
finfo.Variables(3).Dimensions(2).Name      = 'atom';
finfo.Variables(3).Dimensions(2).Length    = natom;
finfo.Variables(3).Dimensions(2).Unlimited = false;
finfo.Variables(3).Dimensions(3).Name      = 'frame';
finfo.Variables(3).Dimensions(3).Length    = nframe;
finfo.Variables(3).Dimensions(3).Unlimited = true;
finfo.Variables(3).Size                    = [3 natom nframe];
finfo.Variables(3).Datatype                = 'single';
finfo.Variables(3).Attributes.Name         = 'units';
finfo.Variables(3).Attributes.Value        = 'angstrom';
finfo.Variables(3).ChunkSize               = [];
finfo.Variables(3).FillValue               = [];
finfo.Variables(3).DeflateLevel            = [];
finfo.Variables(3).Shuffle                 = false;

if exist('box', 'var') && ~isempty(box)
  finfo.Variables(4).Name                 = 'cell_spatial';
  finfo.Variables(4).Dimensions.Name      = 'cell_spatial';
  finfo.Variables(4).Dimensions.Length    = 3;
  finfo.Variables(4).Dimensions.Unlimited = 0;
  finfo.Variables(4).Size                 = 3;
  finfo.Variables(4).Datatype             = 'char';
  finfo.Variables(4).Attributes           = [];
  finfo.Variables(4).ChunkSize            = [];
  finfo.Variables(4).FillValue            = [];
  finfo.Variables(4).DeflateLevel         = [];
  finfo.Variables(4).Shuffle              = false;

  finfo.Variables(5).Name                    = 'cell_angular';
  finfo.Variables(5).Dimensions(1).Name      = 'label';
  finfo.Variables(5).Dimensions(1).Length    = 5;
  finfo.Variables(5).Dimensions(1).Unlimited = false;
  finfo.Variables(5).Dimensions(2).Name      = 'cell_angular';
  finfo.Variables(5).Dimensions(2).Length    = 3;
  finfo.Variables(5).Dimensions(2).Unlimited = false;
  finfo.Variables(5).Size                    = [5 3];
  finfo.Variables(5).Datatype                = 'char';
  finfo.Variables(5).Attributes              = [];
  finfo.Variables(5).ChunkSize               = [];
  finfo.Variables(5).FillValue               = [];
  finfo.Variables(5).DeflateLevel            = [];
  finfo.Variables(5).Shuffle                 = false;

  finfo.Variables(6).Name                    = 'cell_lengths';
  finfo.Variables(6).Dimensions(1).Name      = 'cell_spatial';
  finfo.Variables(6).Dimensions(1).Length    = 3;
  finfo.Variables(6).Dimensions(1).Unlimited = false;
  finfo.Variables(6).Dimensions(2).Name      = 'frame';
  finfo.Variables(6).Dimensions(2).Length    = nframe;
  finfo.Variables(6).Dimensions(2).Unlimited = true;
  finfo.Variables(6).Size                    = [3 nframe];
  finfo.Variables(6).Datatype                = 'double';
  finfo.Variables(6).Attributes.Name         = 'units';
  finfo.Variables(6).Attributes.Value        = 'angstrom';
  finfo.Variables(6).ChunkSize               = [];
  finfo.Variables(6).FillValue               = [];
  finfo.Variables(6).DeflateLevel            = [];
  finfo.Variables(6).Shuffle                 = false;

  finfo.Variables(7).Name                    = 'cell_angles';
  finfo.Variables(7).Dimensions(1).Name      = 'cell_angular';
  finfo.Variables(7).Dimensions(1).Length    = 3;
  finfo.Variables(7).Dimensions(1).Unlimited = false;
  finfo.Variables(7).Dimensions(2).Name      = 'frame';
  finfo.Variables(7).Dimensions(2).Length    = nframe;
  finfo.Variables(7).Dimensions(2).Unlimited = true;
  finfo.Variables(7).Size                    = [3 nframe];
  finfo.Variables(7).Datatype                = 'double';
  finfo.Variables(7).Attributes.Name         = 'units';
  finfo.Variables(7).Attributes.Value        = 'degree';
  finfo.Variables(7).ChunkSize               = [];
  finfo.Variables(7).FillValue               = [];
  finfo.Variables(7).DeflateLevel            = [];
  finfo.Variables(7).Shuffle                 = false;
end

ncwriteschema(filename, finfo);

%% write data
ncwrite(filename, 'spatial', ['xyz']);
ncwrite(filename, 'time', single([1:nframe]));
ncwrite(filename, 'coordinates', single(trj));
if exist('box', 'var') && ~isempty(box)
  ncwrite(filename, 'cell_spatial', ['abc']);
  ncwrite(filename, 'cell_angular', ['alpha'; 'beta '; 'gamma']');
  ncwrite(filename, 'cell_lengths', single(box'));
  ncwrite(filename, 'cell_angles', single(repmat(90, 3, nframe)));
end

