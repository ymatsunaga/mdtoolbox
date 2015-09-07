function make(target)

%% mex files list
mex_opt = cell(1);

% least-squares fitting routine
mex_opt{1} = {'superimpose.c', '-largeArrayDims'};

% mbar auxiliary function;
mex_opt{end+1} = {'mbar_log_wi_jn.c', '-largeArrayDims'};

% 2-d kernel density estimator;
mex_opt{end+1} = {'ksdensity2d.c'};

% 3-d kernel density estimator;
mex_opt{end+1} = {'ksdensity3d.c'};

% reversible transition matrix estimator;
%mex_opt{end+1} = {'msmtransitionmatrix.c', '-largeArrayDims'};

% L-BFGS-B optimization
if ispc
  mex_opt{end+1} = {'lbfgsb_wrapper.c', '-largeArrayDims', '-UDEBUG', '-Ilbfgsb/', ...
                  'lbfgsb/lbfgsb.c','lbfgsb/linesearch.c', ...
                  'lbfgsb/subalgorithms.c','lbfgsb/print.c', ...
                  'lbfgsb/linpack.c','lbfgsb/miniCBLAS.c','lbfgsb/timer.c'};
else
  mex_opt{end+1} = {'lbfgsb_wrapper.c', '-largeArrayDims', '-lm', '-UDEBUG', '-Ilbfgsb/', ...
                  'lbfgsb/lbfgsb.c','lbfgsb/linesearch.c', ...
                  'lbfgsb/subalgorithms.c','lbfgsb/print.c', ...
                  'lbfgsb/linpack.c','lbfgsb/miniCBLAS.c','lbfgsb/timer.c'};
end

nmex = numel(mex_opt);

%% get fullpath to mdtoolbox
fullpath = mfilename('fullpath');
path_mdtoolbox = fileparts(fullpath);

%% change directory
current_dir = pwd;
source_dir = [path_mdtoolbox filesep 'mdtoolbox'];
cd(source_dir);
cleaner = onCleanup(@() cd(current_dir));

if ~isoctave()
  %%%%% MATLAB
  %% default keyword
  if (~exist('target', 'var')) | isempty(target)
    target = 'all';
  end

  %% set options and compile
  if strncmpi(target, 'all', numel('all'))
    for i =1:nmex
      mex_opt{i} = [mex_opt{i}, '-O'];
    end
    compile_mex_matlab(mex_opt);

  elseif strncmpi(target, 'verbose', numel('verbose'))
    for i =1:nmex
      mex_opt{i} = [mex_opt{i}, {'-v', '-O3'}];
    end
    compile_mex_matlab(mex_opt);

  elseif strncmpi(target, 'debug', numel('debug'))
    for i =1:nmex
      mex_opt{i} = [mex_opt{i}, {'-g', '-DDEBUG'}];
    end
    compile_mex_matlab(mex_opt);
    
  elseif strncmpi(target, 'openmp', numel('openmp'))
    for i =1:nmex
      mex_opt{i} = [mex_opt{i}, '-O3', 'LDFLAGS="-fopenmp \$LDFLAGS"', 'CFLAGS="-fopenmp \$CFLAGS"'];
    end
    compile_mex_matlab(mex_opt);

  elseif strncmpi(target, 'clean', numel('clean'))
    delete_mex(mex_opt);

  else
    error(sprintf('cannot find target: %s', target));

  end

else 

  %%%%% Octave
  %% setup
  if (~exist('target', 'var'))
    target = 'all';
  end

  if strncmpi(target, 'all', numel('all'))
    for i = 1:numel(mex_list)
      commandline = sprintf('CFLAGS="-O3" CXXFLAGS="-O3" LDFLAGS="-O3" mkoctfile -v --mex %s', mex_list{i});
      system(commandline);
    end

  elseif strncmpi(target, 'verbose', numel('verbose'))
    for i = 1:numel(mex_list)
      commandline = sprintf('CFLAGS="-O3" CXXFLAGS="-O3" LDFLAGS="-O3" mkoctfile -v --mex %s', mex_list{i});
      system(commandline);
    end

  elseif strncmpi(target, 'debug', numel('debug'))
    for i = 1:numel(mex_list)
      commandline = sprintf('CFLAGS="-g -DDEBUG" CXXFLAGS="-g -DDEBUG" LDFLAGS="-g -DDEBUG" mkoctfile -v --mex %s', mex_list{i});
      system(commandline);
    end
    
  elseif strncmpi(target, 'openmp', numel('openmp'))
    for i = 1:numel(mex_list)
      commandline = sprintf('CFLAGS="-O3 -fopenmp" CXXFLAGS="-O3 -fopenmp" LDFLAGS="-O3 -fopenmp" mkoctfile -v --mex %s -lz -lgomp', mex_list{i});
      system(commandline);
    end

  elseif strncmpi(target, 'clean', numel('clean'))
    delete_mex(mex_list);

  else
    error(sprintf('cannot find target: %s', target));

  end

end

function compile_mex_matlab(mex_opt)
for i = 1:numel(mex_opt)
  fprintf('Compiling MEX code: %s ...\n', mex_opt{i}{1});
  mex(mex_opt{i}{:});
  fprintf('\n');
end

function delete_mex(mex_opt)
for i = 1:numel(mex_opt)
  [~, filename] = fileparts(mex_opt{i}{1});
  mex_exe = [filename '.' mexext];
  fprintf('Deleting MEX exe: %s ...\n', mex_exe);
  if (exist(mex_exe, 'file'))
    delete(mex_exe);
  end

  if isoctave()
    mex_obj = [filename '.o'];
    fprintf('Deleting object file: %s ...\n', mex_obj);
    if (exist(mex_obj, 'file'))
      delete(mex_obj);
    end
  end
end

