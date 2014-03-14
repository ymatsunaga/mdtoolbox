function make(target, opts)

%% mex files list
mex_list{1} = 'superimpose.c'; % least-squares fitting routine
mex_list{2} = 'mbar_log_wi_jn.c'; % mbar auxiliary function;

%% setup
if (nargin < 1) || isempty(target)
  target = 'all';
end

if (nargin < 2)
  opts = '';
end

%% get fullpath to mdtoolbox
fullpath = mfilename('fullpath');
[path_mdtoolbox, filename] = fileparts(fullpath);

%% compile mex files
current_dir = pwd;
source_dir = [path_mdtoolbox filesep 'mdtoolbox'];

cd(source_dir);

if target == 'all'
  opts = '-v'
  compile_mex(mex_list, opts);
elseif target == 'debug'
  option = '-g -DDEBUG';
  compile_mex(mex_list, opts);
elseif target == 'openmp'
  option = 'LDFLAGS="-fopenmp \$LDFLAGS" CFLAGS="-fopenmp \$CFLAGS"';
  compile_mex(mex_list, opts);
elseif target == 'clean'

else
  error(sprintf('cannot find target: %s', target));
end

cd(current_dir);

function compile_mex(mex_list, opts)
for i = 1:numel(mex_list)
  fprintf('Compiling MEX functinon %s ...\n', mex_list{i});
  eval(sprintf('mex %s %s', opts, mex_list{i}));
end

