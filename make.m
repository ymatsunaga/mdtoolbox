function make(target)

%% mex files list
mex_list{1} = 'superimpose.c'; % least-squares fitting routine
mex_list{2} = 'mbar_log_wi_jn.c'; % mbar auxiliary function;
mex_list{3} = 'ksdensity2d.c'; % 2-d kernel density estimator;
mex_list{4} = 'ksdensity3d.c'; % 3-d kernel density estimator;

if ~isoctave()
  %%%%% MATLAB
  %% setup
  if (~exist('target', 'var')) | isempty(target)
    target = 'all';
  end

  %% get fullpath to mdtoolbox
  fullpath = mfilename('fullpath');
  path_mdtoolbox = fileparts(fullpath);

  %% compile mex files
  current_dir = pwd;
  source_dir = [path_mdtoolbox filesep 'mdtoolbox'];

  cd(source_dir);

  if strncmpi(target, 'all', numel('all'))
    opts = '';
    compile_mex(mex_list, opts);

  elseif strncmpi(target, 'verbose', numel('verbose'))
    opts = '-v';
    compile_mex(mex_list, opts);

  elseif strncmpi(target, 'debug', numel('debug'))
    opts = '-g -DDEBUG';
    compile_mex(mex_list, opts);
    
  elseif strncmpi(target, 'openmp', numel('openmp'))
    opts = 'LDFLAGS="-fopenmp \$LDFLAGS" CFLAGS="-fopenmp \$CFLAGS"';
    compile_mex(mex_list, opts);

  elseif strncmpi(target, 'clean', numel('clean'))
    delete_mex(mex_list);

  else
    error(sprintf('cannot find target: %s', target));

  end

  cd(current_dir);
  
else 

  %%%%% Octave
  %% setup
  if (~exist('target', 'var'))
    target = 'all';
  end

  %% get fullpath to mdtoolbox
  fullpath = mfilename('fullpath');
  path_mdtoolbox = fileparts(fullpath);

  %% compile mex files
  current_dir = pwd;
  source_dir = [path_mdtoolbox filesep 'mdtoolbox'];

  cd(source_dir);

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
      commandline = sprintf('CFLAGS="-g" CXXFLAGS="-g" LDFLAGS="-g" mkoctfile -v --mex %s', mex_list{i});
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

  cd(current_dir);

end

function compile_mex(mex_list, opts)
for i = 1:numel(mex_list)
  fprintf('Compiling MEX code: %s ...\n', mex_list{i});
  eval(sprintf('mex %s %s', opts, mex_list{i}));
  fprintf('\n');
end

function delete_mex(mex_list)
for i = 1:numel(mex_list)
  [~, filename] = fileparts(mex_list{i});
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

