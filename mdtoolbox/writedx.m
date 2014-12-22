function writedx(filename, densityMatrix, xi, yi, zi)
%% writedx
% write DX (OpenDX) format file
%
%% Syntax
%# writedx(filename, densityMatrix);
%# writedx(filename, densityMatrix, xi, yi, zi);
% 
%% Description
% This function outputs a DX file for viewing in VMD or PyMOL. It uses
% four inputs to function properly (3D-density data, and grid vectors). While only
% one is necessary (3D-density data), the others position the 3D
% matrix in space. 
%
% * filename      = output dx filename
% * densityMatrix = rectangular-prism 3D matrix holding density values.
%                   All XYZ values correspond to spatial location of
%                   data. [nx x ny x nz double]
% * xi            = grid vector for the x-axis. if omitted, 1:nx is assumed.
%                   [nx x 1 double]
% * yi            = grid vector for the y-axis. if omitted, 1:ny is assumed.
%                   [ny x 1 double]
% * zi            = grid vector for the z-axis. if omitted, 1:nz is assumed.
%                   [nz x 1 double]
% 
%% Example
%# xi = -1:0.1:1; yi = -1:0.1:1; zi = -1:0.1:1;
%# [x, y, z] = meshgrid(xi, yi, zi);
%# densityMatrix = ((x.^2 + y.^2 + z.^2).^0.5);
%# isosurface(x, y, z, densityMatrix, 1.0);
%# writedx('exp.dx', densityMatrix, xi, yi, zi);
%
%% References
% Format of DX (OpenDX) file
% http://www.poissonboltzmann.org/file-formats/mesh-and-data-formats/opendx-scalar-data
% 
%% License
% This routine is based on mat2dx.m originally written by Evan, 
% see README for license information.
% 

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% prepare variables for output function
fprintf('outputting 3D data to %s\n', filename);

voxelLength = zeros(3, 1);
if nargin < 3
  minX = 0;
  minY = 0;
  minZ = 0;
  voxelLength(:) = 1.0;
else
  minX = min(xi);
  minY = min(yi);
  minZ = min(zi);
  voxelLength(1) = xi(2) - xi(1);
  voxelLength(2) = yi(2) - yi(1);
  voxelLength(3) = zi(2) - zi(1);
end

% convenient variables
dimen = size(densityMatrix);
totalElements = prod(dimen);
overFlowVals = mod(totalElements,3);
numRows = floor(totalElements / 3);

% reshape data for fast output
densityMatrix = reshape(permute(densityMatrix, [3 2 1]), [1 totalElements]);
out3DMat = reshape(densityMatrix(1:end - overFlowVals), [3 numRows])';

%% open file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% output data into file
% output header information
fprintf(fid, 'object 1 class gridpositions counts %5d %5d %5d\n', ...
    dimen(1), dimen(2), dimen(3));
fprintf(fid, 'origin %16g %16g %16g\n', minX, minY, minZ);
fprintf(fid, 'delta %16g 0 0\n', voxelLength(1));
fprintf(fid, 'delta 0 %16g 0\n', voxelLength(2));
fprintf(fid, 'delta 0 0 %16g\n', voxelLength(3));
fprintf(fid, 'object 2 class gridconnections counts %5d %5d %5d\n', ...
    dimen(1), dimen(2), dimen(3));
fprintf(fid, 'object 3 class array type double rank 0 items %15d follows\n', ...
    totalElements);

% output density information
newLine = 0;
for n = 1:numRows
    fprintf(fid,'%16E %16E %16E\n',out3DMat(n,:));
    % output status --> often not needed! Function is too fast!
    if ~mod(n, floor(numRows/20))
        fprintf('   %6.0f%%',n/numRows*100);
        newLine = newLine + 1;
        if ~mod(newLine,5)
            fprintf('\n');
        end
    end
end
if overFlowVals > 0 % values not in complete row
    for n = 1:overFlowVals
        fprintf(fid, '%16E',densityMatrix(end-n+1));
    end
    fprintf(fid, '\n');
end
fprintf('\n');

% output file tail
fprintf(fid, ['attribute "dep" string "positions"\n' ...
                'object "regular positions regular connections" class field\n' ...
                'component "positions" value 1\n' ...
                'component "connections" value 2\n' ...
                'component "data" value 3\n']);

%fprintf(fid, 'object "Untitled" call field\n'); % alternate tail

