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
%           step: [nframex1 double]
%           time: [nframex1 double]
%        totener: [nframex1 double]
%          totke: [nframex1 double]
%         energy: [nframex1 double]
%    temperature: [nframex1 double]
%           grms: [nframex1 double]
%        hfctote: [nframex1 double]
%          hfcke: [nframex1 double]
%         ehfcor: [nframex1 double]
%          virke: [nframex1 double]
%          bonds: [nframex1 double]
%         angles: [nframex1 double]
%         urey_b: [nframex1 double]
%      dihedrals: [nframex1 double]
%      impropers: [nframex1 double]
%          cmaps: [nframex1 double]
%        vdwaals: [nframex1 double]
%           elec: [nframex1 double]
%         hbonds: [nframex1 double]
%            asp: [nframex1 double]
%           user: [nframex1 double]
%           vire: [nframex1 double]
%           viri: [nframex1 double]
%         presse: [nframex1 double]
%         pressi: [nframex1 double]
%         volume: [nframex1 double]
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
is_INFO_found = false;
label = {};
ilabel = 0;

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if strncmp(line, 'INFO:', numel('INFO:'))
    [~, line] = strtok(line);
    while ~isempty(line)
      [tmp, line] = strtok(line);
      tmp = strrep(tmp, '-', '_'); % '-' cannnot be used as a filed name in MATLAB
      ilabel = ilabel + 1;
      label{ilabel} = lower(tmp);
    end
    is_INFO_found = true;
    break
  end
end

%% parse data
data = [];
data_remd1 = [];
data_remd2 = [];

fseek(fid, 0, 'bof');
if is_INFO_found

  while ~feof(fid)
    line = strtrim(fgetl(fid));

    if strncmp(line, 'INFO:', numel('INFO:'))
      if ~strncmp(line(13:end), 'STEP', numel('STEP'))
        [~, line] = strtok(line);
        ene_line = cellfun(@str2num, regexp(strtrim(line), '\s*', 'split'));
        data = [data; ene_line];
      end
    end
  end

else
    
  while ~feof(fid)
    line = strtrim(fgetl(fid));

    if strncmp(line, 'RepIDtoParmID:', numel('RepIDtoParmID:'))
      [~, line] = strtok(line);
      data_line = cellfun(@str2num, regexp(strtrim(line), '\s*', 'split'));
      data_remd1 = [data_remd1; data_line];
    end

    if strncmp(line, 'ParmIDtoRepID:', numel('ParmIDtoRepID:'))
      [~, line] = strtok(line);
      data_line = cellfun(@str2num, regexp(strtrim(line), '\s*', 'split'));
      data_remd2 = [data_remd2; data_line];
    end
  end

end

ene = struct;

if ~isempty(data)
  % delete zero-th frame
  %data(1, :) = [];
  for ilabel = 1:numel(label)
    ene = setfield(ene, label{ilabel}, data(:, ilabel));
  end
end

if ~isempty(data_remd1)
  ene = setfield(ene, 'repid2parmid', data_remd1);
end

if ~isempty(data_remd2)
  ene = setfield(ene, 'parmid2repid', data_remd2);
end

