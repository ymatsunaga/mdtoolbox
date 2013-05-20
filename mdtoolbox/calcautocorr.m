function [acf, lag] = calcautocorr(x)
%% calcautocorr
% calculate the sample autocorrelation function of given univariate time-series
%
%% Syntax
%# 
%
%% Description
% If the given input is a matrix, 
% acf is calculated for each column.
%
% * x     - time-series data 
%           [nstep x n double]
% * acf   - sample autocorrelation
%           [nstep x n double]
% * lag   - lag time corresponding to acf
%           [nstep x 1 double]
%
%% Example
%#
%
%% See alo
% 
% 

%% setup
[nstep, nvar] = size(x);

%% sample autocorrelation function
% decenter
x_mean = mean(x);
x = bsxfun(@minus, x, x_mean);

% normalize
x_std = std(x, 1);
x = bsxfun(@rdivide, x, x_std);

% autocorrelation using FFT
lag = 0:floor((nstep-1)/2);
acf = [];
for ivar = 1:nvar
  % FFT using the length of power of 2
  nfft = 2^nextpow2(2*nstep - 1);
  x_k = fft(x, nfft);
  % Wiener–Khinchin theorem theorem
  acf_each = ifft(x_k .* conj(x_k))./nstep;
  % keep values of lag = 0:floor(N/2)
  acf = [acf acf_each(lag+1)];
end


