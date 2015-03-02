function trj_interpolated = interpolatetrj(trj, nsnapshot)
%% interpolatetrj
% interpolate input trajectory by cubic spline
%
%% Syntax
%# trj_interpolated = interpolatetrj(trj, nsnapshot);
%
%% Description
% This routine interpolates infrequently saved trajectory
% by using natural cubic splines.
% This interpolation is just for nice visual inspection. 
% thus should not be treated for scientific analysis. 
%
% * trj              - trajectory data [nstep x natom3]
% * nsnapshot        - the number of snapshots created by the interpolation [scalar integer]
% * trj_interpolated - interpolated trajectory data [nsnapshot x natom3]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# trj_interpolated = interpolatetrj(trj, 1000);
%
%% See alo
% 

%% initialization
if ~exist('nsnapshot', 'var')
  fprintf('interpolated by 100 snapshots\n');
  nsnapshot = 100;
end
nstep = size(trj, 1);
index_time = 1:nstep;
index_time_interpolation = linspace(1, nstep, nsnapshot);

%% natural cubic spline interpolation
cs = spline(index_time, trj');
trj_interpolated = ppval(cs, index_time_interpolation);
trj_interpolated = trj_interpolated';

