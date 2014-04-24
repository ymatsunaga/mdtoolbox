function trj_interpolated = interpolatetrj(trj, nsnapshot)
%% interpolatetrj
% interpolate input trajectory by the natural cubic splines
%
%% Syntax
%# trj_interpolated = interpolatetrj(trj, nsnapshot);
%
%% Description
% This routine interpolates infrequently saved trajectory
% by using the cubic splines.
% This interpolation is just for convenience in visual inspection, 
% and should be carefull for the interpreation of the
% interplated trajectory.
%
% * trj       - coordinates of atoms [nstep x natom3]
% * nsnapshot - the number of snapshots created by the interpolation [scalar integer]
% * trj       - interpolated trajectory [nsnapshot x natom3]
%
%% Example
%# trj = readnetcdf('ala.nc');
%# trj = trj(1:10, :);
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

