function r = isoctave()
%% isoctave
% function that checks if we are in octave
%

persistent x;

if (isempty (x))
  x = exist('OCTAVE_VERSION', 'builtin');
end

r = x;

