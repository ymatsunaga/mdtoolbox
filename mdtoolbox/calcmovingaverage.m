function x_ave = calcmovingaverage(x, nwindow)
%% calcmovingaverage
% calc moving average of input data
% 
%% Syntax
%# x_ave = calcmovingaverage(x);
%# x_ave = calcmovingaverage(x, nwindow);
%
%% Description
% This code calculate a moving average of input data x.
% Average is performed over the rows. 
% This code is expected to be used for a quick check of 
% the trends in the data. 
%
% * x         - input data
%               [nstep x 1 double] or [nstep x m double]
% * nwindow   - size of average window in steps. Default is 10. 
%               [scalar integer]
% 
%% Example
%# e = readamberout('run.out');
%# eptot_ave = calcmovingaverage(e.eptot);
%# plot([e.eptot eptot_ave]);
%# xlabel('step'); ylabel('potential energy [kcal/mol]')
%
%% See also
% 
% 

%% setup
nstep = size(x, 1);
nvar  = size(x, 2);

if ~exist('nwindow', 'var') || isempty(nwindow)
  nwindow = 10
end

%% calculation
x_ave = zeros(nstep, nvar);
a = 1;
b = ones(nwindow, 1)./nwindow;

for istep = 1:nstep
  istart = istep - nwindow;
  if istart < 1
    istart = 1;
  end
  iend = istep + nwindow;
  if iend > nstep
    iend = nstep;
  end
  index = istart:iend;
  x_ave(istep, :) = mean(x(index, :));
end

% x_ave = filter(b, a, x);
% for ivar = 1:nvar
%   %x_ave(:, ivar) = conv(x(:, ivar), mask, 'same');
%   x_ave(:, ivar) = filter(b, a, x(:, ivar));
% end

