function [index, dmin] = assignvoronoi(ref, data)
%% assignvoronoi
% assign Voronoi cell indices to the elements of data by checking the distances between ref and data
%
%% Syntax
%# index = assignvoronoi(ref, data)
%# [index, dmin] = assignvoronoi(ref, data)
%
%% Description
% * ref      - generating points of Voronoi cells [n x ndim double]
% * data     - some trajectory or data set to be assigned [nstep x ndim double]
% * index    - indices of Voronoi cell. indices correspond to the rows of ref [nstep double]
% * box      - distance between the data to the nearest generating point [nstep double]
%
%% Example
%# ref  = randn(100, 1);
%# data = randn(10000, 1);
%# [index, dmin] = assignvoronoi(ref, data);
% 
%% See also
% assign1dbins, assign2dbins, 
%

%% setup
[ncell, ndim1] = size(ref);
[nstep, ndim2] = size(data);
assert(ndim1 == ndim2, 'dimensions of ref and data do not match...');

%% calculation
index = zeros(nstep, 1);
dmin = zeros(nstep, 1);

for istep = 1:nstep
  dev = bsxfun(@minus, ref, data(istep, :));
  dist = sum(dev.^2,2);
  [dmin_each, index_each] = min(dist);
  index(istep) = index_each;
  dmin(istep) = sqrt(dmin_each);
end


