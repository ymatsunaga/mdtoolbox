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
regexp_getlabel = '[^(1-9\(\)\s-=\.#)]*';
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
      match_token = regexp(line, regexp_gettoken, 'match');
      for i = 1:numel(match_token)
        token = match_token{i};
        match_label = regexp(token, regexp_getlabel, 'match');
        ilabel = ilabel + 1;
        label{ilabel} = lower(match_label{1});
      end
      line = strtrim(fgetl(fid));
    end

    break;
  
  end % end of if line == NSTEP

end
nlabel = ilabel;


%% parse data
frewind(fid);
ilabel = 1;
data = [];
data_each = zeros(1, nlabel);

while ~feof(fid)
  line = strtrim(fgetl(fid));

  if strfind(line, 'NMR restraints:')
    continue;
  end

  if strfind(line, 'A V E R A G E S   O V E R')
    break;
  end

  if strfind(line, 'R M S  F L U C T U A T I O N S')
    break;
  end

  % parse one step data
  if strncmp(line, 'NSTEP', numel('NSTEP'))

    while ~(strncmp(line, '---', numel('---')))

      match_token = regexp(line, regexp_gettoken, 'match');

      for i = 1:numel(match_token)

        token = match_token{i};
        match_label = regexp(token, regexp_getlabel, 'match');
        match_data = regexp(token, regexp_getdata, 'match');

        for ilabel = 1:nlabel
          if strncmpi(match_label{1}, label{ilabel}, numel(label{ilabel}))
            data_each(ilabel) = str2num(match_data{1}(2:end));
          end
        end
        
      end

      line = strtrim(fgetl(fid));
    
    end
  
  end

  % store one step data
  if ilabel == nlabel
    data = [data; data_each];
    ilabel = 1;
    data_each = zeros(1, nlabel);
  end
  
end

for i = 1:nlabel
  if strcmp(label{i}, 'nb')
    label{i} = 'nb14';
  elseif strcmp(label{i}, 'eel')
    label{i} = 'eel14';
  end
end

ene = struct;
for ilabel = 1:nlabel
  ene = setfield(ene, label{ilabel}, data(:, ilabel));
end

