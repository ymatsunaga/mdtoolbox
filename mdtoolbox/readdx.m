function [densityMatrix, xi, yi, zi] = readdx(filename)
%% readdx
% read DX (OpenDX) format file
%
%% Syntax
%# densityMatrix = readdx(filename);
%# [densityMatrix, xi, yi, zi] = readdx(filename);
%
%% Description
% This function reads in a DX file and translates it into a 3D matrix in
% matlab. DX files are file formats for storing 3D data in plaintext, and
% they are readable by PyMOL and VMD. One inputs the DX file path, and the
% program outputs a structure with all the information about the DX
% 
% * filename      = input dx filename
% * densityMatrix = rectangular-prism 3D matrix holding density values.
%                   All XYZ values correspond to spatial location of
%                   data. [nx x ny x nz double]
% * xi            = grid vector for the x-axis
%                   [nx x 1 double]
% * yi            = grid vector for the y-axis
%                   [ny x 1 double]
% * zi            = grid vector for the z-axis
%                   [nz x 1 double]
%
%% Example
%# [densityMatrix, xi, yi, zi] = readdx('water_density.dx');
%# [x, y, z] = meshgrid(xi, yi, zi);
%# isosurface(x, y, z, densityMatrix, 3.0);
%
%% References
% Format of DX (OpenDX) file
% http://www.poissonboltzmann.org/file-formats/mesh-and-data-formats/opendx-scalar-data
% 
%% License
% This routine is based on dx2mat.m originally written by Evan, 
% see README for license information.
% 

%% read every line in file
assert(ischar(filename), 'Please specify valid filename as the first argument')
fprintf('reading in file %s\n', filename);
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% read all into memory
rawText = fread(fid, inf, '*char');

%% interpret introductory data
fprintf('importing data about 3D matrix\n');

% parse lines by end-of-lines
splitLines = strread(rawText, '%s', 'delimiter', '\n');
numLines = length(splitLines);

% use introductory data to deduce file length and grid positions
voxelLength = [0 0 0]';
lineCounter = 1;
while ( lineCounter < numLines )
    
    % find number of voxels in each dimension
    if ~isempty(strmatch('object 1 class gridpositions counts', splitLines(lineCounter)))
        
        % decode string
        strResults = regexp(splitLines(lineCounter), '(\d+)','match');
        numResults = str2double(strResults{1}(:));
        
        % save voxel counts
        numXvoxels = numResults(2);
        numYvoxels = numResults(3);
        numZvoxels = numResults(4);
    end
    
    % find origin or data
    if ~isempty(strmatch('origin', splitLines(lineCounter)))
        
        % decode string
        %strResults = regexp(splitLines(lineCounter), '(\d+)','match');
        strResults = regexp(splitLines(lineCounter), '([-+]?\d+(\.\d+)?)','match');
        numResults = str2double(strResults{1}(:));
        
        % save voxel counts
        minX = numResults(1);
        minY = numResults(2);
        minZ = numResults(3);
    end
    
    % find voxel length in each dimension
    if ~isempty(strmatch('delta ', splitLines(lineCounter)))
        
        % decode string
        %strResults = regexp(splitLines(lineCounter), '(\d+)','match');
        strResults = regexp(splitLines(lineCounter), '([-+]?\d+(\.\d+)?)','match');
        numResults = str2double(strResults{1}(:));
        
        % record voxel dimensions
        voxelLength = voxelLength + numResults;
    end
    
    % find total number of data points
    if ~isempty(strmatch('object 3 class array type double rank 0 items', splitLines(lineCounter)))
        
        % decode string, save results
        strResults = regexp(splitLines(lineCounter), '(\d+)','match');
        totalDXelements = str2double(strResults{1}(3));
        
        firstLineOfData = lineCounter + 1;
    end
    
    % exit while loop
    if (exist('numXvoxels', 'var') && sum(voxelLength ~= 0) == 3 && ...
            exist('minX', 'var') && exist('totalDXelements', 'var'))
        lineCounter = numLines;
    end
    
    lineCounter = lineCounter + 1;
    
end

% check number of 3D data elements
if ( (numXvoxels * numYvoxels * numZvoxels) ~= totalDXelements )
    fprintf('ERROR!\n\tdimensional mismatch in number of 3D elements!\n');
    return;
end

%% interpret 3D grid data
fprintf('reading in 3D data\n');

% pull in main bulk of data
linearMatrix = zeros(1, totalDXelements);
lineCounter = firstLineOfData;
lastLineOfData = firstLineOfData + ceil(totalDXelements / 3) - 1;
matCounter = 1;
newLine = 0;
while ( lineCounter < lastLineOfData )
    
    strResults = regexp( splitLines(lineCounter), '(\S+)','match');
    numResults = str2double(strResults{1}(:));
    
    linearMatrix(matCounter:matCounter+2) = numResults;
    
    lineCounter = lineCounter + 1;
    matCounter  = matCounter + 3;
    
    if ~mod(lineCounter, floor(totalDXelements/20))
        fprintf('   %6.0f%%',lineCounter/ceil(totalDXelements / 3)*100);
        newLine = newLine + 1;
        if ~mod(newLine,5)
            fprintf('\n');
        end
    end
    
end

% pull in last line of data
strResults = regexp( splitLines(lastLineOfData), '(\S+)','match');
numResults = str2double(strResults{1}(:));
linearMatrix(matCounter:end) = numResults;

fprintf('      100%%\n');

%% prepare output
fprintf('finishing ...\n');

% save information
xi = minX + ((1:numXvoxels) - 1)*voxelLength(1);
yi = minY + ((1:numYvoxels) - 1)*voxelLength(2);
zi = minZ + ((1:numZvoxels) - 1)*voxelLength(3);

% 3D matrix
densityMatrix = permute(reshape(linearMatrix, [numZvoxels numYvoxels numXvoxels]), [3 2 1]);

