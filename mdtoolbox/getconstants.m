function constant = getconstants()

constant.KB  = 0.00198719168260038; % Boltzmann constant (kB)in kcal/mol/K
constant.T   = 300.0;               % temperature in K
constant.KBT = constant.KB * constant.T; % kBxT in kcal/mol

%% TODO transformation of velocity unit

