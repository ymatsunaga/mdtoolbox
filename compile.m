%% This script compiles MEX files

% mex files list
mex_list{1} = 'superimpose.c'; % least-squares fitting routine

% get fullpaths to mdtoolbox
this_file = mfilename('fullpath');
[path_mdtoolbox, filename, ext] = fileparts(this_file);

% compile mex files
current_dir = pwd;
source_dir = [path_mdtoolbox filesep 'mdtoolbox'];

cd(source_dir);

for i = 1:numel(mex_list)
  fprintf('Compiling MEX functinon %s ...\n', mex_list{i});
  mex(mex_list{i});
end

cd(current_dir);

% configure MATLAB path for mdtoolbox
addpath([path_mdtoolbox filesep 'mdtoolbox']);
%addpath([path_mdtoolbox filesep 'test']);

% cleaning
clear mex_list;
clear this_file;
clear path_mdtoolbox filename, ext;
clear current_dir source_dir;

