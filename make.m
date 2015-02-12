function make(target)

%% mex files list
mex_list{1} = 'superimpose.c'; % least-squares fitting routine
mex_list{2} = 'mbar_log_wi_jn.c'; % mbar auxiliary function;
mex_list{3} = 'ksdensity3d.c'; % 3-d radial distribution function;

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
end

