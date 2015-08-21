function [c2, trj_index2] = msmtarjan(c, trj_index)
%% msmtarjan
% Tarjan's algorithm
%
%% Syntax
%# [c2, index2] = msmtarjan(c, index);
%

nstate = size(c, 1);
if nstate > 1000
  nRecursionLimit = get(0,'RecursionLimit');
  set(0,'RecursionLimit', nstate*2);
end

index_counter = 1;
stack = [];
for i = 1:nstate
  lowlinks{i} = [];
  index{i} = [];
end
result = [];

for i = 1:nstate
  graph{i} = find(c(i, :));
end

for i = 1:nstate
  if isempty(lowlinks{i})
    [index_counter, stack, lowlinks, index, result] = strongconnect(graph, i, index_counter, stack, lowlinks, index, result);
  end
end
c2 = result;
[~, id_ergodic] = max(cellfun(@(x) numel(x), c2));
index_ergodic = sort(c2{id_ergodic});
nstate_new = numel(index_ergodic);
c2 = c(index_ergodic, index_ergodic);

if nargin > 1
  if iscell(trj_index)
    ncell = numel(trj_index);
    for icell = 1:ncell
      id = ismember(trj_index{icell}, index_ergodic);
      trj_index{icell}(~id) = NaN;
      for istate = 1:nstate_new
        id = (trj_index{icell} == index_ergodic(istate));
        trj_index{icell}(id) = istate;
      end
    end
  else
    id = ismember(trj_index, index_ergodic);
    trj_index(~id) = NaN;
    for istate = 1:nstate_new
      id = (trj_index == index_ergodic(istate));
      trj_index(id) = istate;
    end
  end
  trj_index2 = trj_index;
else
  trj_index2 = [];
end

if nstate > 1000
  set(0,'RecursionLimit', nRecursionLimit);
end


%% this resursion is based on http://www.logarithmic.net/pfh-files/blog/01208083168/tarjan.py
function [index_counter, stack, lowlinks, index, result] = strongconnect(graph, node, index_counter, stack, lowlinks, index, result)
%
%

index{node} = index_counter;
lowlinks{node} = index_counter;
index_counter = index_counter + 1;
stack(end+1) = node;

successors = graph{node};

for successor = successors
  if isempty(lowlinks{successor})
    [index_counter, stack, lowlinks, index, result] = strongconnect(graph, successor, index_counter, stack, lowlinks, index, result);
    lowlinks{node} = min(lowlinks{node}, lowlinks{successor});
  elseif ismember(successor, stack)
    lowlinks{node} = min(lowlinks{node}, index{successor});
  end
end

if lowlinks{node} == index{node}
  connected_component = [];
  
  while true
    successor = stack(end);
    stack(end) = [];
    connected_component(end+1) = successor;
    if successor == node
      break;
    end
  end
  
  result{end+1} = connected_component;
end

