function ene = readamberout(filename)
%% readamberout
% read amber output file
%
%% Syntax
%# ene = readamberout(filename);
%
%% Description
% Output(ene) is a structure variable. 
%
% * fname     - filename of amber output
% * ene       - structure data 
%        nstep: [nstepx1 double]
%         time: [nstepx1 double]
%         temp: [nstepx1 double]
%        press: [nstepx1 double]
%         etot: [nstepx1 double]
%        ektot: [nstepx1 double]
%        eptot: [nstepx1 double]
%         bond: [nstepx1 double]
%        angle: [nstepx1 double]
%        dihed: [nstepx1 double]
%         nb14: [nstepx1 double]
%        eel14: [nstepx1 double]
%      vdwaals: [nstepx1 double]
%        eelec: [nstepx1 double]
%       ehbond: [nstepx1 double]
%    restraint: [nstepx1 double]
%       eamber: [nstepx1 double]
% 
%% Example
%# ene = readamberout('run.log');
%# plot(ene.eptot)
%
%% See also
% readnamdout, readgenesisout
% 

%% define regular expression text
regexp_gettoken = '\w[^=]*=\s*[-\d\.]+';
regexp_getlabel = '[^(0-9\(\)\s-=\.)]*';
regexp_getdata  = '=\s*[-s\d\.]+';

%% open file
fid = fopen(filename, 'r');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

%% parse label
label = {};
ilabel = 0;

while 1
  line = strtrim(fgetl(fid));

  if strncmp(line, 'NSTEP', numel('NSTEP'))

    while ~(strncmp(line, '---', numel('---')))
      match = regexp(line, regexp_gettoken, 'match');
      for i = 1:numel(match)
        token = match{i};
        match2 = regexp(token, regexp_getlabel, 'match');
        ilabel = ilabel + 1;
        label{ilabel} = lower(match2{1});
      end
      line = strtrim(fgetl(fid));
    end

    break;
  
  end % enf of if line == NSTEP

end

for i = 1:numel(label)
  if strcmp(label{i}, 'nb')
    label{i} = 'nb14';
  elseif strcmp(label{i}, 'eel')
    label{i} = 'eel14';
  end
end

%% parse data
frewind(fid);
ilabel = 0;
data = [];
data_each = zeros(1, numel(label));

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if regexp(line, 'A V E R A G E S   O V E R')
    break;
  end

  if regexp(line, 'R M S  F L U C T U A T I O N S')
    break;
  end

  isdata = false;
  if regexp(line, regexp_gettoken)
    for i = 1:numel(label)
      if regexpi(line, ['^' label{i} '[\s\(]'])
        isdata = true;
      end
    end
    if strncmpi(line, '1-4 ', numel('1-4 '))
      isdata = true;
    end
  end

  if isdata
    match = regexp(line, regexp_gettoken, 'match');
    for i = 1:numel(match)
      token = match{i};
      match2 = regexp(token, regexp_getdata, 'match');
      ilabel = ilabel + 1;
      data_each(ilabel) = str2num(match2{1}(2:end));
    end
  end

  if ilabel == numel(label)
    data = [data; data_each];
    ilabel = 0;
  end
  
end

ene = struct;
for ilabel = 1:numel(label)
  ene = setfield(ene, label{ilabel}, data(:, ilabel));
end

