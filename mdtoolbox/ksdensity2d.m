function [f, xi, yi] = ksdensity2d(data, xi, yi, weight, bw)
%% [f, xi, yi] = ksdensity2d(data, xi, yi, weight, bw)
%
%

%% setup
nstep = size(data, 1);

if ~exist('xi', 'var') || isempty(xi)
  xi = linspace(min(data(:, 1)), max(data(:, 1)), 100);
elseif iscolumn(xi)
  xi = xi';
end
nx = numel(xi);

if ~exist('yi', 'var') || isempty(yi)
  yi = linspace(min(data(:, 2)), max(data(:, 2)), 100);
elseif iscolumn(yi)
  yi = yi';
end
ny = numel(yi);

if ~exist('weight', 'var') || isempty(weight)
  weight = ones(nstep, 1);
end
weight = weight./sum(weight);

if ~exist('bw', 'var') || isempty(bw)
  bw = zeros(2, 1);
  
  sig = median(abs(xi - median(xi)))/0.6745;
  if sig <= 0, sig = max(xi) - min(xi); end
  if sig > 0
    bw(1) = sig * (1/nstep)^(1/6);
  else
    bw(1) = 1;
  end

  sig = median(abs(yi - median(yi))) / 0.6745;
  if sig <= 0, sig = max(yi) - min(yi); end
  if sig > 0
    bw(2) = sig * (1/nstep)^(1/6);
  else
    bw(2) = 1;
  end
end
bw

%% compute the kernel density estimates
diffx2 = (bsxfun(@minus, data(:, 1), xi)./bw(1)).^2;
diffy2 = (bsxfun(@minus, data(:, 2), yi)./bw(2)).^2;
gaussx = exp(-0.5 * diffx2)./(sqrt(2*pi).*bw(1));
gaussy = exp(-0.5 * diffy2)./(sqrt(2*pi).*bw(2));

f = zeros(ny, nx);
for istep = 1:nstep
  t2 = bsxfun(@times, gaussx(istep, :), gaussy(istep, :)');
  f = f + t2*weight(istep);
end

