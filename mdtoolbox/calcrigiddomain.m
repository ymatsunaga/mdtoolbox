function [indexOfDomain, F_history] = calcrigiddomain(ref, pmode, variance, ndomain, alpha)
%% calcrigiddomain
% calculate rigid-domains from mode vectors
%
%% Syntax
%# indexOfDomain = calcrigiddomain(ref, pmode, variance);
%# indexOfDomain = calcrigiddomain(ref, pmode, variance, ndomain);
%# indexOfDomain = calcrigiddomain(ref, pmode, variance, ndomain, alpha);
%
%% Description
% This routine performs the Quasi-Rigid Domain Decomposition by
% using the method proposed by Micheletti and coworkers
% (Biophys. J. 2009) from given input mode vectors and the
% variances in those directions. 
% Simulated annealing algorithm is used for the minimization of a
% objective function F.
%
%% Example
%# trj = readdcd('ak_ca.dcd');
%# [ref, trj] = meanstructure(trj);
%# [p, pmode, variances] = calcpca(trj);
%# [indexOfDomain, F_history] = calcrigiddomain(ref, pmode(:, 1:20), variances(1:20), 5);
% In PyMOL,
% PyMOL> @domain.pml
%
%% See also
% anm calcpca
% 
%% References
% [1] R. Potestio, F. Pontiggia, and C. Micheletti, Biophys. J. 96, 4993 (2009).
% [2] T. Aleksiev, R. Potestio, F. Pontiggia, S. Cozzini, and
%     C. Micheletti, Bioinformatics 25, 2743 (2009). 
% 

%% constants
natom = size(ref, 2)/3;
nmode = size(pmode, 2);
Rc = 7.0;
temperature = 100000;

if ~exist('ndomain', 'var') || isempty(ndomain)
  ndomain = 5;
end

if ~exist('alpha', 'var') || isempty(alpha)
  alpha = 10.0;
end

%% elasticity
elasticity = zeros(natom, natom);
for i = 1:natom
  for j = 1:natom
    index_i = (3*(i-1)+1):(3*i);
    index_j = (3*(j-1)+1):(3*j);
    for imode = 1:nmode
      delta_mode = pmode(index_i, imode) - pmode(index_j, imode);
      delta_distance = ref(index_i) - ref(index_j);
      elasticity(i, j) = elasticity(i, j) + variance(imode)*(dot(delta_mode, delta_distance).^2);
    end
  end
end

%% connectivity
connectivity = calcdistancematrix(ref);
connectivity = 0.5 * (1.0 + tanh(Rc - connectivity));
for iatom = 1:natom
  connectivity(iatom, iatom) = 0.0;
end

%% objective function F
indexOfDomain = randi(ndomain, natom, 1);
F = calcF(indexOfDomain, elasticity, connectivity, ndomain, alpha, Rc);

%% simulated annealing (minimize F)
F_best = inf;

icount = 1;
F_history = F;
iatom = 1;
while true
  % trial move
  iatom = randomwalk(iatom, natom);
  idomain = indexOfDomain(iatom);

  %jatom = randomwalk(iatom, natom);
  %idomain_try = indexOfDomain(jatom);

  idomain_try = randi(ndomain);

  %idomain_try = mod(idomain, ndomain) + 1;

  % elasticity
  delta_F = 0.0;
  index     = find(indexOfDomain == idomain);
  index_try = find(indexOfDomain == idomain_try);
  delta_F = delta_F - 0.5*sum(elasticity(iatom, index));
  delta_F = delta_F + 0.5*sum(elasticity(iatom, index_try));

  % connectivity
  delta_F = delta_F + 0.5*sum(connectivity(iatom, index));
  delta_F = delta_F - 0.5*sum(connectivity(iatom, index_try));

  % judgement
  if exp(-delta_F/temperature) > rand
    F = F + delta_F;
    indexOfDomain(iatom) = idomain_try;
  end

  % preprocess
  if mod(icount, 10000) == 0
    fprintf('Trial = %d  T = %f  F = %f\n', icount, temperature, F);
  end
  temperature = 0.99999*temperature; %exponential scheduling
  icount = icount + 1;
  F_history = [F_history; F];
  if F < F_best
    F_best = F;
    indexOfDomain_best = indexOfDomain;
  end
  
  if temperature < 0.1
    break;
  end
end
indexOfDomain = indexOfDomain_best;

%% output pymol script
%F = calcF(indexOfDomain, elasticity, connectivity, ndomain, alpha, Rc);
fid = fopen('domain.pml', 'w');
c = linspace(0, 1, ndomain);
cmap = colormap;
ncolor = size(cmap, 1);
c = c * (ncolor-1) + 1;
cs = spline(1:ncolor, cmap');
rgb = ppval(cs, c);
rgb = rgb';
rgb(rgb < 0) = 0;
rgb(rgb > 1) = 1;
for idomain = 1:ndomain
  fprintf(fid, 'set_color domain%d, [%f, %f, %f]\n', idomain, rgb(idomain, 1), rgb(idomain, 2), rgb(idomain, 3));
end
fprintf(fid, '\n');
for idomain = 1:ndomain
  index = find(indexOfDomain == idomain);
  if ~isempty(index)
    fprintf(fid, 'select domain%d, resi ', idomain);
    for i = 1:(numel(index)-1)
      fprintf(fid, '%d+', index(i));
    end
    fprintf(fid, '%d\n', index(end));
  end
end
fprintf(fid, '\n');
for idomain = 1:ndomain
  fprintf(fid, 'color domain%d, domain%d\n', idomain, idomain);
end
fprintf(fid, '\n');
fclose(fid);


%% trial move by random walk
function iatom_out = randomwalk(iatom, natom)
% random walk
if rand > 0.5
  iatom_out = iatom + 1;
else
  iatom_out = iatom - 1;
end
% periodic boundary
if iatom_out <= 0
  iatom_out =  iatom_out + natom;
elseif iatom_out > natom
  iatom_out =  iatom_out - natom;
end


%% calculate objective function F
function F = calcF(indexOfDomain, elasticity, connectivity, ndomain, alpha, Rc)
F = 0;
for idomain = 1:ndomain
  index = find(indexOfDomain == idomain);
  if numel(index) >= 2
    pair = nchoosek(index, 2);
    F = F + 0.5*sum(arrayfun(@(x,y) elasticity(x,y), pair(:, 1), pair(:, 2)));
  end
end

for idomain = 1:ndomain
  index_i = find(indexOfDomain == idomain);
  for jdomain = (idomain+1):ndomain
    index_j = find(indexOfDomain == jdomain);
    if (numel(index_i) >= 1) && (numel(index_j) >= 1)
      [p, q] = meshgrid(index_i(:), index_j(:));
      pair = [p(:) q(:)];
      F = F + alpha*0.5*sum(arrayfun(@(x,y) connectivity(x,y), pair(:, 1), pair(:, 2)));
    end
  end
end

