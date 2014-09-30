function ene = readgenesisout(filename)
%% readgenesisout
% read genesis output file 
%
%% Syntax
%# ene = readgenesisout(filename);
%
%% Description
% Output(ene) is a structure variable. 
%
% * fname     - filename of genesis output
% * ene       - structure data 
%           step: [nstepx1 double]
%           time: [nstepx1 double]
%        totener: [nstepx1 double]
%          totke: [nstepx1 double]
%         energy: [nstepx1 double]
%    temperature: [nstepx1 double]
%           grms: [nstepx1 double]
%        hfctote: [nstepx1 double]
%          hfcke: [nstepx1 double]
%         ehfcor: [nstepx1 double]
%          virke: [nstepx1 double]
%          bonds: [nstepx1 double]
%         angles: [nstepx1 double]
%         urey_b: [nstepx1 double]
%      dihedrals: [nstepx1 double]
%      impropers: [nstepx1 double]
%          cmaps: [nstepx1 double]
%        vdwaals: [nstepx1 double]
%           elec: [nstepx1 double]
%         hbonds: [nstepx1 double]
%            asp: [nstepx1 double]
%           user: [nstepx1 double]
%           vire: [nstepx1 double]
%           viri: [nstepx1 double]
%         presse: [nstepx1 double]
%         pressi: [nstepx1 double]
%         volume: [nstepx1 double]
% 
%% Example
%# ene = readgenesisout('run.log');
%# plot(ene.energy)
%
%% See also
% readamberout, readnamdout
% 

%% define regular expression text
regexp_getdata  = '(\S+?\.\S{5}|\d+)';

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse label
label = {};
ilabel = 0;

while 1
  line = strtrim(fgetl(fid));

  if strncmp(line, 'INFO:', numel('INFO:'))
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

  if strncmp(line, 'INFO:', numel('INFO:'))
    [~, line] = strtok(line);
    ene_line = cellfun(@str2num, regexp(strtrim(line), '\s*', 'split'));
    data = [data; ene_line];
  end
end

% delete zero-th step
%data(1, :) = [];

ene = struct;
for ilabel = 1:numel(label)
  ene = setfield(ene, label{ilabel}, data(:, ilabel));
end

