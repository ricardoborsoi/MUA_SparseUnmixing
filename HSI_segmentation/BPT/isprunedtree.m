function b = isprunedtree(T)

% checks if the tree T has already been pruned. If it is the case, it
% returns a true boolean, otherwise false.

b = isfield(T(1),'pruning');

% security check
if ~isequal(b,isfield(T(end),'pruning'))
    error('Field name ''pruning'' not consistent along the tree structure array')
end
