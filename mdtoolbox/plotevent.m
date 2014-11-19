function gobj = plotevent(event_time)
%% plotevent
% plot poisson events by stems
%
%% Syntax
%# plotcolor(event_time)
%
%% Description
% * event_time - XXXXXXXX.
% * gobj       - graphics object (formerly called as 'handle graphics')
%
%% Example
%# [time_event, time_interval] = simulatepoisson(1, 20);
%# plotevent(time_event)
%

%% preparation
nevent = numel(event_time);
yvalue = ones(1, nevent);

%% plot stems
gobj = stem(event_time, yvalue, 'k-', 'Marker', 'none', 'LineWidth', 1.5);
axis([0 max(event_time) 0 2])
formatplot

xlabel('time', 'fontsize', 25);


