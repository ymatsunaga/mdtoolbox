function trj = wrap(trj, box, resname)
%% wrap
% wrap waters and ions into the primary box under Periodic Boundary Conditions (PBCs)
%
%% Syntax
%# trj = wrap(trj, box);        % ignore residue-information
%# trj = wrap(trj, box, residue_name); % consider residue
%
%% Description
%
% * trj          - trajectory [nframe x natom3 double]
% * box          - box size [nframe x 3 double]
% * residue_name - residue name lists [nframe x n array chars]
%
%% Example
%# [trj, box] = readnetcdf('run.nc');
%# trj2 = wrap(trj, box);
%# writenetcdf('run_wrap', trj2, box);
%
%% See also
%
%% References
%

%% setup
nframe = size(trj, 1);

%% setup cell

if ~exist('resname', 'var')
  for iframe = 1:nframe
    trj(iframe, 1:3:end) = trj(iframe, 1:3:end) + box(iframe, 1)*0.5;
    trj(iframe, 2:3:end) = trj(iframe, 2:3:end) + box(iframe, 2)*0.5;
    trj(iframe, 3:3:end) = trj(iframe, 3:3:end) + box(iframe, 3)*0.5;
    trj(iframe, 1:3:end) = trj(iframe, 1:3:end) - floor(trj(iframe, 1:3:end)./box(iframe, 1))*box(iframe, 1);
    trj(iframe, 2:3:end) = trj(iframe, 2:3:end) - floor(trj(iframe, 2:3:end)./box(iframe, 2))*box(iframe, 2);
    trj(iframe, 3:3:end) = trj(iframe, 3:3:end) - floor(trj(iframe, 3:3:end)./box(iframe, 3))*box(iframe, 3);
    trj(iframe, 1:3:end) = trj(iframe, 1:3:end) - box(iframe, 1)*0.5;
    trj(iframe, 2:3:end) = trj(iframe, 2:3:end) - box(iframe, 2)*0.5;
    trj(iframe, 3:3:end) = trj(iframe, 3:3:end) - box(iframe, 3)*0.5;
  end
else
  % wrap waters
  index = find(selectname(resname, 'WAT') | selectname(resname, 'TIP3'));
  trj_sub = trj(:, to3(index));
  for iframe = 1:nframe
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) + box(iframe, 1)*0.5;
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) + box(iframe, 2)*0.5;
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) + box(iframe, 3)*0.5;
    x_cell = floor(trj_sub(iframe, 1:9:end)./box(iframe, 1));
    y_cell = floor(trj_sub(iframe, 2:9:end)./box(iframe, 2));
    z_cell = floor(trj_sub(iframe, 3:9:end)./box(iframe, 3));
    x_cell = kron(x_cell, [1 1 1]);
    y_cell = kron(y_cell, [1 1 1]);
    z_cell = kron(z_cell, [1 1 1]);
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) - x_cell*box(iframe, 1);
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) - y_cell*box(iframe, 2);
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) - z_cell*box(iframe, 3);
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) - box(iframe, 1)*0.5;
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) - box(iframe, 2)*0.5;
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) - box(iframe, 3)*0.5;
  end
  trj(:, to3(index)) = trj_sub;
  % wrap ions
  index = find(selectname(resname, 'SOD') | selectname(resname, 'CL') ...
             | selectname(resname, 'CLA') | selectname(resname, 'MG') ...
             | selectname(resname, 'POT') | selectname(resname, 'CES') ...
             | selectname(resname, 'CAL') | selectname(resname, 'ZN2') ...
             | selectname(resname, 'Na+') | selectname(resname, 'Cl-'));
  trj_sub = trj(:, to3(index));
  for iframe = 1:nframe
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) + box(iframe, 1)*0.5;
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) + box(iframe, 2)*0.5;
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) + box(iframe, 3)*0.5;
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) - floor(trj_sub(iframe, 1:3:end)./box(iframe, 1))*box(iframe, 1);
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) - floor(trj_sub(iframe, 2:3:end)./box(iframe, 2))*box(iframe, 2);
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) - floor(trj_sub(iframe, 3:3:end)./box(iframe, 3))*box(iframe, 3);
    trj_sub(iframe, 1:3:end) = trj_sub(iframe, 1:3:end) - box(iframe, 1)*0.5;
    trj_sub(iframe, 2:3:end) = trj_sub(iframe, 2:3:end) - box(iframe, 2)*0.5;
    trj_sub(iframe, 3:3:end) = trj_sub(iframe, 3:3:end) - box(iframe, 3)*0.5;
  end
  trj(:, to3(index)) = trj_sub;
end

