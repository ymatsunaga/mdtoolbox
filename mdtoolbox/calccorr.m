function [corr, corr_atom] = calccorr(trj, lagtime)
%% calccovar
% calculate (time-lagged) cross-correlation matrix from input trajectory
%
%% Syntax
%# corr = calccorr(trj);
%# corr = calccorr(trj, lagtime);
%# [corr, corr_atom] = calccorr(trj, lagtime)
%
%% Description
% This routine calculates the Pearson correlation coefficients 
% between atomi and atom j, and make a cross-correlation matrix
% from input trajectory
%
% * trj       - trajectory of coordinates [nframe x 3natom]
% * lagtime   - lag time in the unit of frames. The default is 0.
%               [scalar integer]
% * corr      - time-lagged cross correlation matrix of anisotropic fluctuations. 
%               [double natom x natom]
% * corr_atom - time-lagged cross correlation matrix of atomic fluctuations. 
%               [double natom x natom]
%
%% Example
%# trj = readdcd('ak.dcd');
%# [~, trj] = meanstructure(trj);
%# [~, corr_atom] = calccorr(trj);
%# imagesc(corr_atom, [-1, 1]);
%# axis xy; axis square; colorbar;
%# xlabel('residue', 'fontsize', 25);
%# ylabel('residue', 'fontsize', 25);
%# exportas('corr');
%
%% See also
% calccovar
% 

%% setup
%nframe = size(trj, 1);
natom3 = size(trj, 2);
natom = natom3/3;

if ~exist('lagtime', 'var') || isempty(lagtime)
  lagtime = 0;
end

%% calc covariance matrix
if nargout == 1
  covar = calccovar(trj, lagtime);
  if lagtime ~= 0
    var0 = var(trj);
  end
elseif nargout > 1
  [covar, covar_atom] = calccovar(trj, lagtime);
  if lagtime ~= 0
    var0 = var(trj);
    var0_atom = zeros(1, natom);
    for i = 1:natom
      var0_atom(i) = sum(var0((3*(i-1)+1):(3*i)));
    end
  end
end

%% normalize covariance matrix
corr = zeros(natom3, natom3);
for i = 1:natom3
  for j = 1:natom3
    if lagtime == 0
      corr(i,j) = covar(i,j)./sqrt(covar(i,i)*covar(j,j));
    else
      corr(i,j) = covar(i,j)./sqrt(var0(i)*var0(j));
    end
  end
end

if nargout > 1
  corr_atom = zeros(natom, natom);
  for i = 1:natom
    for j = 1:natom
      if lagtime == 0
        corr_atom(i,j) = covar_atom(i,j)./sqrt(covar_atom(i,i)*covar_atom(j,j));
      else
        corr_atom(i,j) = covar_atom(i,j)./sqrt(var0_atom(i)*var0_atom(j));
      end
    end
  end
end

