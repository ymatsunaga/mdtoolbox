function p_i = msmensemble(emission_im, h_m, p0_i)
%% msmensemble
% this code estimate equilibrium probabilities p_i from projected distribution data 
%
%% Syntax
%# p_i = msmensemble(emission_im, h_m);
%# p_i = msmensemble(emission_im, h_m, p0_i);
%
%% Description
% This code estimate equilibrium probabilities p_i from projected
% distribution data  
%
% * emission_im - probability of observing grid m given state i
% * h_m         - projected distribution data
% * p0_i        - initial estimate for equilibrium probabilities
% * p_i         - estimated equilibrium probabilities
%
%% Example
%# [state, observation] = msmgenerate(10^4, T, emission_im, p_i);
%# h_m = zeros(1, size(emission_im, 2));
%# for m = 1:size(emission_im, 2)
%#   h_m(m) = sum(observation == m);
%# end
%# p_i = msmensemble(emission_im, h_m, p0_i)
%
%
%% See also
% msmbaumwelch, msmbaumwelchdb
%

%% iterative algorithm
nstate = size(emission_im, 1);
ngrid = size(emission_im, 2);

posterior_im = zeros(nstate, ngrid);
if exist('index', 'var') && ~isempty(index)
  p_i = p0_i;
else
  p_i = ones(1, nstate)./nstate;
end
p_i_old = Inf(1, nstate);

while max(abs(p_i_old - p_i)) > 10^(-6)
  p_i_old = p_i;
  for i = 1:nstate
    posterior_im(i, :) = p_i(i) * emission_im(i, :) ./ sum(bsxfun(@times, p_i', emission_im));
  end
  for i = 1:ncluster
    index = ~isnan(posterior_im(i, :));
    p_i(i) = sum(h_m{itry}(index) .* posterior_im(i, index)) ./ sum(h_m{itry});
  end
end

