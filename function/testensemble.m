function [ratio, ratio_fit, center, a_fit, a_theory, a_boot] = testensemble(e1, e2, t1, t2)
%% testensemble
% test ensemble
% 
%% Example
%
%% References
% [1] M. R. Shirts, J. Chem. Theory Comput. 9, 909 (2013).
% 

%% setup
nbin = 100;
N1 = numel(e1);
N2 = numel(e2);

min_e1 = min(e1);
min_e2 = min(e2);
min_e  = max([min_e1 min_e2]);

max_e1 = max(e1);
max_e2 = max(e2);
max_e  = min([max_e1 max_e2]);

s = getconstants();
a_theory(1) = 0;
a_theory(2) = 1./(s.KB*t2) - 1./(s.KB*t1);

%% probability ratio
edge = linspace(min_e, max_e + nbin*eps, nbin + 1);
center = edge + 0.5*(edge(2)-edge(1));
center(end) = [];

c1 = histc(e1, edge);
c1(end) = [];
c1 = c1./sum(c1);

c2 = histc(e2, edge);
c2(end) = [];
c2 = c2./sum(c1);

ratio = log(c2./c1);

%% maximum likelihood
a_fit = maximum_likelihood_estimation(e1, e2);
ratio_fit = - a_fit(1) - a_fit(2)*center;

%% plot
hold off
plot(center, ratio, 'b-');
hold on;
plot(center, ratio_fit, 'r-');

%% bootstrap
nboot = 100;
a_boot = zeros(nboot, 2);

for i = 1:nboot
  e1_boot = e1(randi(numel(e1), numel(e1), 1));
  e2_boot = e2(randi(numel(e2), numel(e2), 1));
  a_tmp = maximum_likelihood_estimation(e1_boot, e2_boot);
  a_boot(i, :) = a_tmp;
end

%% print the result
fprintf('Theory: beta2 - beta1 = %f\n', a_theory(2));
fprintf('Sim:    beta2 - beta1 = %f +- %f\n', a_fit(2), std(a_boot(:,2)));

function a = maximum_likelihood_estimation(e1, e2)
%a = fminsearch(@(a) negative_likelihood(a, e1, e2), [0 0], optimset('MaxIter', 10^10, 'TolFun', 10^(-10)));
a = fminsearch(@(a) negative_likelihood(a, e1, e2), [0 0]);

function y = negative_likelihood(a, e1, e2)
N1 = numel(e1);
N2 = numel(e2);
M = log(N1/N2);
%y = (sum(log(f(- M - a(1) - a(2)*e1))) + sum(log(f(+ M + a(1) + a(2)*e2))))./(N1+N2);
y = (sum(log(f(- M - a(1) - a(2)*e1))) + sum(log(f(+ M + a(1) + a(2)*e2))));

function y = f(x)
y = 1 + exp(x);

