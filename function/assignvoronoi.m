function [bin_n, dmin_n] = assignvoronoi(ref, data_n)
%% assignvoronoi
%
%% Syntax
%#
%
%% Description
%
%% Example
%#
% 
%% See also
%
%% References
%

[ncell,ndim1] = size(ref);
[nstep,ndim2] = size(data_n);

if ndim1 ~= ndim2
  error('ndim1 and ndim2 is not equal...');
end
ndim = ndim1;

bin_n = zeros(nstep,1);
dmin_n = zeros(nstep,1);
for istep = 1:nstep
  tdiff = bsxfun(@minus,ref,data_n(istep,:));
  d = sum(tdiff.^2,2);
  [dmin,bin] = min(d);
  bin_n(istep) = bin;
  dmin_n(istep) = sqrt(dmin);
end

% for i = 1:nstep
%   d_min = -1.0;
%   icell = 0;
%   for j = 1:ncell
%     d = sum((data_n(i,:) - ref(j,:)).^2);
%     if (d < d_min) || (icell == 0)
%       d_min = d;
%       icell = j;
%     end
%   end
%   bin_n(i) = icell;
%   dmin_n(i) = sqrt(d_min);
% end


