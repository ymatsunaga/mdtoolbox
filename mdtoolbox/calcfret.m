function [E, NA, ND] = calcfret(distance, T, N, R0)
%% calcfret
% calculate FRET efficiency and donor and acceptor photons from one-dimensional distance time-series
%
%% Syntax
%# [E, NA, ND] = calcfret(distance, T, N, R0)
%
%% Description
%
% * distance   - distance between donor and acceptor dyes (in Angstrom)
% * T          - length of time-bin (in frames)
% * N          - total photon counting rate per time-bin
% * R0         - Forster radius (in Angstrom)
% * E          - FRET efficiencies in time-bins
% * NA         - Acceptor photon countings in time-bins
% * ND         - Donor photon countings in time-bins
%
%% Example
%# T = 1000;
%# N = 10;
%# R0 = 40;
%# nframe = T*100;
%# distance = abs(cumsum(randn(nframe,1)) + R0);
%# [E, NA, ND] = calcfret(distance, T, N, R0);
%# plot(distance)
%# plot(E)
%# 
%
%% See alo
% 
% 

%% setup
nframe = numel(distance);

%% calculate FRET efficiency
iframe = 1;
icount = 1;

while (iframe+T-1) <= nframe
  % generate the total photon counts (N_A + N_B) according to the Poisson distribution
  NA_ND = poissrnd(N);
  % calculate the average FRET efficiency
  e = mean(1./(1 + (distance(iframe:(iframe+T-1))./R0).^6));
  % generate the acceptor photon counts
  NA(icount) = binornd(NA_ND, e);
  ND(icount) = NA_ND - NA(icount);
  E(icount) = NA(icount)./(NA_ND);
  
  % update local variables
  iframe = iframe + T;
  icount = icount + 1;
end

