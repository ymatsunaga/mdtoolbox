function [p_ci_best, values] = clusterinformation(x, temperature, nc, nReplicates)
%% clusteringbyinformation
% clusterize samples according to an information-based criterion
%
%% Syntax
%# [p_ci, s] = clusteringbyinformation(x, temperature, numCluster, nInitialCondition);
%
%% Description
% assuming p(i) = 1/N;
%
%% Example
%# 
%# [p_ci, s] = clusteringbyinformation(q, 1/30, 3, 10);
%# scatter(q(:, 1), q(:, 2), 50, p_ci, 'fill')
%
%% See alao
% clusteringbykmeans, clusteringbykcenter
%
%% References
% Slonim, N., Atwal, G. S., Tkačik, G. & Bialek, W. 
% "Information-based clustering." PNAS 102, 18297–18302 (2005).
%
%% TODO
% support for p(i) ~= 1/N
%

%% constants
maxIteration = 2000;
tolerance = 10^(-10);

%% setup
[nframe, ndim] = size(x);

if ~exist('temperature', 'var')
    temperature = 1/30;
end

if ~exist('nc', 'var')
    nc = round(sqrt(ndim));
end

if ~exist('nReplicates', 'var')
    nReplicates = 10;
end

%% clustering
values.f           = -inf;
values.similarity  = 0;
values.compression = 0;

p_ci_best = zeros(nframe, nc);
p_ci_old  = zeros(nframe, nc);
work      = zeros(nframe, nc);

for ireplica = 1:nReplicates
  % initialize
  p_ci = rand(nframe, nc);
  p_ci = bsxfun(@rdivide, p_ci, sum(p_ci, 2));
  p_c  = sum(p_ci) ./ nframe;
  
  % calc similarities s(c) and s(C; i)
  s_c  = zeros(1, nc);
  s_ci = zeros(nframe, nc);
  for i = 1:nframe
    similarity_i = calcsimilarity(x(i,:), x);
    s_ci(i, :)   = similarity_i * p_ci;
  end
  % end
  s_c  = sum(p_ci .* s_ci);
  s_ci = bsxfun(@rdivide, s_ci, p_c) ./ nframe;
  s_c  = s_c .* (1./nframe.^2) .* (1./p_c.^2);
  
  % solve self-consistent equation
  f           = -inf;
  delta       = inf;
  similarity  = 0;
  compression = 0;
  icount      = 1;

  while delta > tolerance
    % update p_ci
    p_ci_old = p_ci;
    p_ci = bsxfun(@minus, 2*s_ci, s_c) ./ temperature;
    p_ci = bsxfun(@minus, p_ci, max(p_ci, [], 2)); % scaling to avoid overflow
    p_ci = exp(p_ci);
    p_ci = bsxfun(@times, p_c, p_ci);
    p_ci = bsxfun(@rdivide, p_ci, sum(p_ci, 2));
    p_c  = sum(p_ci) ./ nframe;

    % calc similarity <s>
    for i = 1:nframe
      similarity_i = calcsimilarity(x(i,:), x);
      s_ci(i, :)   = similarity_i * p_ci;
    end
    s_c  = sum(p_ci .* s_ci);
    s_ci = bsxfun(@rdivide, s_ci, p_c) ./ nframe;
    s_c  = s_c .* (1./nframe.^2) .* (1./p_c.^2);
    similarity = sum(p_c .* s_c);

    % calc compression I
    work = bsxfun(@rdivide, p_ci, p_c);
    work(work <= eps) = 1.0;
    compression = sum(sum( p_ci .* log2(work) ));
    compression = compression ./ nframe;

    % calc free energy f
    f = similarity - temperature * compression;

    delta = max(max(abs(p_ci_old - p_ci)));
    
    icount = icount + 1;
    if icount > maxIteration
      disp(sprintf('maximum number (= %d) of iteration is reached', maxIteration));
      break
    end
  end

  if f > values.f
    values.f           = f;
    values.similarity  = similarity;
    values.compression = compression;
    p_ci_best          = p_ci;
  end
  
  disp(sprintf('%d iteration, f = %f', ireplica, values.f));
end


% function ss = calcsimilarity(xi, x)
%   [iframe, ndim] = size(xi);
%   [nframe, ndim] = size(x);
%   ss = zeros(iframe, nframe);
%   cc = cell(ndim, 1);
%   for idim = 1:ndim
%     cc{idim} = bsxfun(@minus, xi(:, idim), x(:, idim)');
%   end
%   for idim = 1:ndim
%     ss = ss + cc{idim}.^2;
%   end
%   ss = -sqrt(ss);


% function ss = calcsimilarity(xi, x)
%   [iframe, ndim] = size(xi);
%   [nframe, ndim] = size(x);
%   ss = zeros(iframe, nframe);
%   cc = zeros(iframe, nframe, ndim);
%   for idim = 1:ndim
%     cc(:, :, idim) = bsxfun(@minus, xi(:, idim), x(:, idim)');
%   end
%   ss = -sqrt(sum(cc.^2, 3));


% function ss = calcsimilarity(xi, x)
%   [iframe, ndim] = size(xi);
%   [nframe, ndim] = size(x);
%   a = reshape(xi, iframe, 1, ndim);
%   b = reshape(x, 1, nframe, ndim);
%   ss = - sqrt(sum((a(:, ones(nframe, 1), :) - b(ones(iframe, 1), :, :)).^2, 3));


function ss = calcsimilarity(xi, x)
  dx = bsxfun(@minus, x, xi);
  ss = sum(dx.^2, 2);
  ss = -sqrt(ss');


% function ss = calcsimilarity(xi, x)
%   xi = xi - mean(xi);
%   x = bsxfun(@minus, x, mean(x, 2));
%   xi_std = std(xi);
%   x_std = std(x, 0, 2);
%   ss = sum(bsxfun(@times, x, xi), 2);
%   ss = ss ./ x_std;
%   ss = ss ./ xi_std;
%   ss = abs(ss');

