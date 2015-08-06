function rmsf = calcrmsf(trj, index, nblock)
%% calcrmsd
% calculate root mean square fluctuatons (RMSFs) of atomic positions
%
%% Syntax
%# rmsf = calcrmsf(trj);
%# rmsf = calcrmsf(trj, index_atom);
%# rmsf = calcrmsf(trj, index_atom, nblock);
%# rmsf = calcrmsf(trj, [], nblock);
%
%% Description
% Calculate RMSF of atomic positions
%
% * trj        - trajectory [nframe x (natom*3) double]
% * index_atom - atom index or logical index specifying atoms to be fitted to the average structure
% * nblock     - the number of blocks for block averaging
% * rmsf       - root mean square fluctuatons (RMSF)
%                if nblock > 1, twice of standard deviation is added in 2nd column for 95% confidence interval
%
%% Example
%# trj = readdcd('trj_ca.dcd');
%# rmsf = calcrmsf(trj);
%# plot(rmsf); xlabel('residue'); ylabel('RMSF [Angstrom]')
%#
%# rmsf = calcrmsf(trj, [], 5);
%# errorbar(1:nresidue, rmsf(:, 1), rmsf(:, 2)); xlabel('residue'); ylabel('RMSF [Angstrom]')
%
%% See alo
% superimpose, meanstructure
% 

%% preparation
nframe = size(trj, 1);
natom = size(trj, 2)/3;

% index
if ~exist('index', 'var') || isempty(index)
  index = 1:natom;
end

% nblock
if ~exist('nblock', 'var') || isempty(nblock)
  nblock = 1;
end

%% evaluate RMSF
interface = round(linspace(0, nframe, nblock+1));
rmsf = {};
for iblock = 1:nblock
  istart = interface(iblock)+1;
  iend = interface(iblock+1);
  if nblock > 1
    fprintf('[block %d] from frame %d to frame %d\n', iblock, istart, iend);
  end
  rmsf{iblock} = kernelfunction(trj(istart:iend, :), index);
end

rmsf = cell2mat(rmsf);
if nblock > 1
  rmsf = [mean(rmsf, 2), 2*std(rmsf, [], 2)];
end

%%%%%%%% kernel function
function rmsf = kernelfunction(trj, index)
[ref, trj] = meanstructure(trj, index);
trj_var = var(trj);
rmsf = sqrt(trj_var(1:3:end) + trj_var(2:3:end) + trj_var(3:3:end))';

