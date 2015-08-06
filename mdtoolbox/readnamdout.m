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
%           ts: [nframex1 double]
%         bond: [nframex1 double]
%        angle: [nframex1 double]
%        dihed: [nframex1 double]
%        imprp: [nframex1 double]
%        elect: [nframex1 double]
%          vdw: [nframex1 double]
%     boundary: [nframex1 double]
%         misc: [nframex1 double]
%      kinetic: [nframex1 double]
%        total: [nframex1 double]
%         temp: [nframex1 double]
%    potential: [nframex1 double]
%       total3: [nframex1 double]
%      tempavg: [nframex1 double]
%     pressure: [nframex1 double]
%    gpressure: [nframex1 double]
%       volume: [nframex1 double]
%     pressavg: [nframex1 double]
%    gpressavg: [nframex1 double]
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
      tmp = strrep(tmp, '-', '_'); % '-' cannnot be used as a filed name in MATLAB
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

% delete zero-th frame
data(1, :) = [];

ene = struct;
for ilabel = 1:numel(label)
  ene = setfield(ene, label{ilabel}, data(:, ilabel));
end

