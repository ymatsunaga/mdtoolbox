function interarrival_time = simulatepoisson(time_width, rate)
%% simulatepoisson
% simulates inter-arrival times of the Poisson process within the given time-width and rate
%

%% simulation
t_cumsum = 0;
interarrival_time = [];

while true
  t = - (1./rate) * log(rand);
  t_cumsum = t_cumsum + t;
  if t_cumsum > time_width
    break
  end
  interarrival_time = [interarrival_time; t];
end

