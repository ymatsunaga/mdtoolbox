function ene = readnamdout(filename)
%% readnamdout
% read namd output file
%
%% Syntax
%# ene = readnamdout(filename);
%
%% Description
% Output(ene) is a structure variable. 
%
% * fname     - filename of namd output
% * ene       - structure data 
%           ts: [nstepx1 double]
%         bond: [nstepx1 double]
%        angle: [nstepx1 double]
%        dihed: [nstepx1 double]
%        imprp: [nstepx1 double]
%        elect: [nstepx1 double]
%          vdw: [nstepx1 double]
%     boundary: [nstepx1 double]
%         misc: [nstepx1 double]
%      kinetic: [nstepx1 double]
%        total: [nstepx1 double]
%         temp: [nstepx1 double]
%    potential: [nstepx1 double]
%       total3: [nstepx1 double]
%      tempavg: [nstepx1 double]
%     pressure: [nstepx1 double]
%    gpressure: [nstepx1 double]
%       volume: [nstepx1 double]
%     pressavg: [nstepx1 double]
%    gpressavg: [nstepx1 double]
% 
%% Example
%# ene = readnamdout('run.log');
%# plot(ene.potential)
%
%% See also
% readamberout, readgenesisout
% 

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse label
label = {};
ilabel = 0;

while 1
  line = strtrim(fgetl(fid));

  if strncmp(line, 'ETITLE:', numel('ETITLE:'))
    [~, line] = strtok(line);
    while ~isempty(line)
      [tmp, line] = strtok(line);
      tmp = strrep(tmp, '-', '_'); % '-' cannnot be used as a filed name
      ilabel = ilabel + 1;
      label{ilabel} = lower(tmp);
    end
    break
  end
end

%% parse data
data = [];

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if strncmp(line, 'ENERGY:', numel('ENERGY:'))
    [~, line] = strtok(line);
    ene_line = cellfun(@str2num, regexp(strtrim(line), '\s*', 'split'));
    data = [data; ene_line];
  end
end

% delete zero-th step
data(1, :) = [];

ene = struct;
for ilabel = 1:numel(label)
  ene = setfield(ene, label{ilabel}, data(:, ilabel));
end


