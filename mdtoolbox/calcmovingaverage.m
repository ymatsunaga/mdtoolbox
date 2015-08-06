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
%               [nframe x 1 double] or [nframe x m double]
% * nwindow   - size of average window in frames. Default is 10. 
%               [scalar integer]
% 
%% Example
%# e = readamberout('run.out');
%# eptot_ave = calcmovingaverage(e.eptot);
%# plot([e.eptot eptot_ave]);
%# xlabel('frame'); ylabel('potential energy [kcal/mol]')
%
%% See also
% 
% 

%% setup
nframe = size(x, 1);
nvar  = size(x, 2);

if ~exist('nwindow', 'var') || isempty(nwindow)
  nwindow = 10
end

%% calculation
x_ave = zeros(nframe, nvar);
a = 1;
b = ones(nwindow, 1)./nwindow;

for iframe = 1:nframe
  istart = iframe - nwindow;
  if istart < 1
    istart = 1;
  end
  iend = iframe + nwindow;
  if iend > nframe
    iend = nframe;
  end
  index = istart:iend;
  x_ave(iframe, :) = mean(x(index, :));
end

% x_ave = filter(b, a, x);
% for ivar = 1:nvar
%   %x_ave(:, ivar) = conv(x(:, ivar), mask, 'same');
%   x_ave(:, ivar) = filter(b, a, x(:, ivar));
% end

