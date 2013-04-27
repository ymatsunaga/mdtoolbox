%% define umbrella window centers
window_center = 0:3:180;
nwindow = numel(window_center);

%% read dihedral angle data
dihedral = {};
for i = 1:nwindow
  filename = sprintf('../prod/run_%d.dat', window_center(i));
  x = load(filename);
  dihedral{i} = x(:, 2);
end

%% define a function handle of bias energy for each umbrella window
fhandle = {};
for i = 1:nwindow
  k = 0.5*0.12184;
  fhandle{i} = @(x) k*(periodic(x, window_center(i))).^2;
end

s = getconstants;
kbt = s.KB*300;

%% calculate probability along the dihedral angle and calculate
%% the potential of mean force (PMF) along the dihedral angle
[f, prob, center, h, bias, ndata] = wham(linspace(-1, 181, 82), fhandle, dihedral, kbt);
pmf = - kbt * log(prob);
pmf = pmf - pmf(1);

%% plot the PMF
hold off
plot(center, pmf, 'k-');
formatplot
xlabel('angle [degree]', 'fontsize', 20);
ylabel('PMF [kcal/mol]', 'fontsize', 20);
axis([-1 181 -8 12]);
exportas('script_pmf')
hold off

%% save results
save -v7.3 script_pmf.mat;

