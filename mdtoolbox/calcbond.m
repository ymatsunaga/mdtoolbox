function bond = calcbond(trj, pair)
%% calcbond
% calculate distances between pair atom from their Cartesian coordinates
%
%% Syntax
%# bond = calcbond(trj);
%
%% Description
% Calculate distances from the input trajectory of Cartesian coordinates.
% Pairs, whose distances are calculated, can be specified via the
% variable (pair).
%
% * trj    - coordinates of atoms [nstep x natom3]
% * pair   - pair indices whose distances are calculated [npair x 2]
% * bond   - distances between the pairs [nstep x npair]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# pair = [5 6; 9 10];
%# bond = calcbond(trj, pair);
%
%% See alo
% calcangle, calcdihedral
% 

%% initialization
if ~exist('pair', 'var')
  %fprintf('calculating the distances between the 1st and 2nd atoms\n');
  pair = [1 2];
end

nstep = size(trj, 1);
npair = size(pair, 1);

%% calculation
bond = zeros(nstep, npair);

for ipair = 1:npair
  index1 = to3(pair(ipair, 1));
  index2 = to3(pair(ipair, 2));
  bond(:, ipair) = sqrt(sum((trj(:, index1) - trj(:, index2)).^2, 2));
end

