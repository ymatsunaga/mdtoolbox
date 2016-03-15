function [mi, mi_excess] = calcmutinf(x, y, k)
%% calcmutinf
% evaluate mutual information
%
%% Syntax
%# mi = calcmutinf(x, y, k)
%# [mi, mi_excess] = calcmutinf(x, y, k)
%
%% Description
% 
% * x - d-dimensional data [double nframe x d]
% * y - d-dimensional data [double nframe x d]
% * k - k-nearest neighbor
%
% * mi - mutual information (I1 in Ref.) (in unit of natrual logarithm)
% * mi_excess - excess mutual information (in unit of natrual logarithm)
%
%% Example
%# 
%
%% See also
% 
% 
%% References
% This function using the method described in
% A. Kraskov, H. Stögbauer, and P. Grassberger
% Phys. Rev. E 69, 066138 (2004).
%

%% setup
nframe_x = size(x, 1);
nframe_y = size(y, 1);

d_x = size(x, 2);
d_y = size(y, 2);

assert(nframe_x == nframe_y, 'frames of x and y do not mtach...');
assert(d_x == d_y, 'dimensions of x and y do not match...');

nframe = nframe_x;
d = d_x;

if ~exist('k', 'var') || isempty(k)
  k = 1;
end

%% mutual information
mi = kernel_func(x, y, k);

%% excess mutual information
if nargout > 1
  mi_perm = zeros(10, 1);
  for i = 1:10
    y_dummy = y(randperm(nframe), :);
    mi_perm(i) = kernel_func(x, y_dummy, k);
  end
  mi_excess = mi - mean(mi_perm);
else
  mi_excess = [];
end

%% kernel function
function mi = kernel_func(x, y, k)
nframe = size(x, 1);
nx = zeros(nframe, 1);
ny = zeros(nframe, 1);
for iframe = 1:nframe
  dx = sqrt(sum(bsxfun(@minus, x, x(iframe, :)).^2, 2));
  dy = sqrt(sum(bsxfun(@minus, y, y(iframe, :)).^2, 2));
  dz = max(dx, dy);
  dz = sort(dz, 'ascend');
  e = dz(k+1);
  
  nx(iframe) = sum(dx < e);
  ny(iframe) = sum(dy < e);
end
nx = nx - 1;
ny = ny - 1;

mi = psi(k) - sum(psi(nx + 1) + psi(ny + 1))/nframe + psi(nframe);

