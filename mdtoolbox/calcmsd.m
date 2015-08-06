function [msd, time, coefficient] = calcmsd(trj, index, ps_per_frame, nblock);
%% calcmsd
% calculate mean square displacements (MSDs) and diffusion coefficient from trajectory
%
%% Syntax
%# [msd, time, coefficient] = calcrmsf(trj);
%# [msd, time, coefficient] = calcrmsf(trj, index_atom);
%# [msd, time, coefficient] = calcrmsf(trj, index_atom, ps_per_frame);
%# [msd, time, coefficient] = calcrmsf(trj, index_atom, ps_per_frame, nblock);
%
%% Description
% Calculate MSDs and diffusion coefficient from trajectory
%
% * trj          - trajectory [nframe x (natom*3) double]
% * index_atom   - atom index or logical index specifying atoms to be fitted to the average structure
% * ps_per_frame - picosecond per frame
% * nblock       - the number of blocks for block averaging
% * msd          - MSD [nframe2 x 1]
%                  if nblock > 1, twice of standard deviation is added in 2nd column for 95% confidence interval
% * time         - time in ps [nframe2 x 1]
% * msd          - diffusion coefficient in 10^{-5}cm^{2}/s [1 or 2]
%                  if nblock > 1, twice of standard deviation is added in 2nd column for 95% confidence interval
%
%% Example
%# trj = readdcd('run.dcd');
%# psf = readpsf('ionize.psf');
%# index_atom = selectname(psf.atom_name, 'OH2');
%# [msd, time, coefficient] = calcmsd(trj, index_atom, 2);
%# loglog(time, msd);
%
%% See alo
% calcrdf calcrmsf
% 

%% preparation
nframe = size(trj, 1);

% index
natom3 = size(trj, 2);
natom = natom3/3
if ~exist('index', 'var') || isempty(index)
  index = 1:natom;
end
trj = trj(:, to3(index));

% ps_per_frame
if ~exist('ps_per_frame', 'var') || isempty(ps_per_frame)
  if nargout > 2
    error('ps_per_frame is required for calculation of diffusion coefficient');
  else
      ps_per_frame = 2;
  end
end

% nblock
if ~exist('nblock', 'var') || isempty(nblock)
  nblock = 1;
end

%% evaluate MSD and diffusion coefficient
msd = {};
time = [];
coefficient = {};
nframe_of_block = floor(nframe/nblock);
interface = round(linspace(0, nframe_of_block*nblock, nblock+1));
for iblock = 1:nblock
  istart = interface(iblock)+1;
  iend = interface(iblock+1);
  if (nblock > 1)
    fprintf('[block %d] from frame %d to frame %d\n', iblock, istart, iend);
  end
  if (nargout > 2)
    [msd{iblock}, time, coefficient{iblock}] = kernelfunction(trj(istart:iend, :), ps_per_frame);
  else
    [msd{iblock}, time] = kernelfunction(trj(istart:iend, :), ps_per_frame);
  end
end

msd = cell2mat(msd);
if (nblock > 1)
  msd = [mean(msd, 2), 2*std(msd, [], 2)];
end

if (nargout > 2)
  coefficient = cell2mat(coefficient);
  if (nblock > 1)
    coefficient = [mean(coefficient, 2), 2*std(coefficient, [], 2)];
    fprintf('[block-averaging] Diffusion coefficient: D=%4.2f +- %4.2f [10^{-5}cm^{2}/s]\n', coefficient(1), coefficient(2));
  end
end

%%%%%%%% kernel function
function [msd, time, coefficient] = kernelfunction(trj, ps_per_frame)
nframe = size(trj, 1);
natom = size(trj, 2)/3;

nframe_of_msd = floor(nframe/10);
msd = zeros(nframe_of_msd, 1);
time = (1:nframe_of_msd)'*ps_per_frame;
for iframe = 1:nframe_of_msd
  nblock = floor(nframe/iframe);
  d = zeros(1, natom);
  d = d + sum((trj(iframe:iframe:end, 1:3:end) - trj(1:iframe:(end-iframe+1), 1:3:end)).^2, 1);
  d = d + sum((trj(iframe:iframe:end, 2:3:end) - trj(1:iframe:(end-iframe+1), 2:3:end)).^2, 1);
  d = d + sum((trj(iframe:iframe:end, 3:3:end) - trj(1:iframe:(end-iframe+1), 3:3:end)).^2, 1);
  assert(size(trj(iframe:iframe:end,:),1) == nblock, 'nblock is not consistent');
  d = d./nblock;
  msd(iframe) = mean(d, 2);
end

if (nargout > 2)
  % time scale faster than 10 ps is ignored for diffusion coefficient calculation
  index_time = time > 10;
  if(any(~index_time))
    disp('message: time scale faster than 10 ps is ignored for diffusion coefficient calculation.');
  end
  p = polyfit(time(index_time), msd(index_time), 1);
  coefficient = p(1)*10./6;
  fprintf('Diffusion coefficient: D=%4.2f [10^{-5}cm^{2}/s]\n', coefficient);
end

