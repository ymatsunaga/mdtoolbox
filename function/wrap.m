function trj = wrap(trj, box, resname)
%% wrap
% wrap waters and ions into the primary box under Periodic Boundary Conditions (PBCs)
%
%% Syntax
%# trj = wrap(trj, box);        % ignore residue-information
%# trj = wrap(trj, box, resid); % consider residue
%
%% Description
%
% * trj        - trajectory [nstep x natom3 double]
% * box        - box size [nstep x 3 double]
% * resname    - residue name lists [nstep x n array chars]
% * trj2       - trajectory wrapped in to the primary box [nstep x natom3 double]
%
%% Example
%# [trj, box] = readnetcdf('run.dcd');
%# trj2 = wrap(trj, box);
%# writenetcdf('run_wrap', trj2, box);
%
%% See also
%
%% References
%

%% setup
[nstep, natom3] = size(trj);
natom = natom3 / 3;

%% setup cell

if nargin < 3
  for istep = 1:nstep
    trj(istep, 1:3:end) = trj(istep, 1:3:end) + box(istep, 1)*0.5;
    trj(istep, 2:3:end) = trj(istep, 2:3:end) + box(istep, 2)*0.5;
    trj(istep, 3:3:end) = trj(istep, 3:3:end) + box(istep, 3)*0.5;
    trj(istep, 1:3:end) = trj(istep, 1:3:end) - floor(trj(istep, 1:3:end)./box(istep, 1))*box(istep, 1);
    trj(istep, 2:3:end) = trj(istep, 2:3:end) - floor(trj(istep, 2:3:end)./box(istep, 2))*box(istep, 2);
    trj(istep, 3:3:end) = trj(istep, 3:3:end) - floor(trj(istep, 3:3:end)./box(istep, 3))*box(istep, 3);
    trj(istep, 1:3:end) = trj(istep, 1:3:end) - box(istep, 1)*0.5;
    trj(istep, 2:3:end) = trj(istep, 2:3:end) - box(istep, 2)*0.5;
    trj(istep, 3:3:end) = trj(istep, 3:3:end) - box(istep, 3)*0.5;
  end
else
  % wrap waters
  index = find(selectname(resname, 'WAT') | selectname(resname, 'TIP3'));
  trj_sub = trj(:, to3(index));
  for istep = 1:nstep
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) + box(istep, 1)*0.5;
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) + box(istep, 2)*0.5;
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) + box(istep, 3)*0.5;
    x_cell = floor(trj_sub(istep, 1:9:end)./box(istep, 1));
    y_cell = floor(trj_sub(istep, 2:9:end)./box(istep, 2));
    z_cell = floor(trj_sub(istep, 3:9:end)./box(istep, 3));
    x_cell = kron(x_cell, [1 1 1]);
    y_cell = kron(y_cell, [1 1 1]);
    z_cell = kron(z_cell, [1 1 1]);
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) - x_cell*box(istep, 1);
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) - y_cell*box(istep, 2);
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) - z_cell*box(istep, 3);
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) - box(istep, 1)*0.5;
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) - box(istep, 2)*0.5;
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) - box(istep, 3)*0.5;
  end
  trj(:, to3(index)) = trj_sub;
  % wrap ions
  index = find(selectname(resname, 'SOD') | selectname(resname, 'CL') ...
             | selectname(resname, 'CLA') | selectname(resname, 'MG') ...
             | selectname(resname, 'POT') | selectname(resname, 'CES') ...
             | selectname(resname, 'CAL')  | selectname(resname, 'ZN2'));
  trj_sub = trj(:, to3(index));
  for istep = 1:nstep
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) + box(istep, 1)*0.5;
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) + box(istep, 2)*0.5;
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) + box(istep, 3)*0.5;
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) - floor(trj_sub(istep, 1:3:end)./box(istep, 1))*box(istep, 1);
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) - floor(trj_sub(istep, 2:3:end)./box(istep, 2))*box(istep, 2);
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) - floor(trj_sub(istep, 3:3:end)./box(istep, 3))*box(istep, 3);
    trj_sub(istep, 1:3:end) = trj_sub(istep, 1:3:end) - box(istep, 1)*0.5;
    trj_sub(istep, 2:3:end) = trj_sub(istep, 2:3:end) - box(istep, 2)*0.5;
    trj_sub(istep, 3:3:end) = trj_sub(istep, 3:3:end) - box(istep, 3)*0.5;
  end
  trj(:, to3(index)) = trj_sub;
end


