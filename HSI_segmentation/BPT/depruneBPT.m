function tree = depruneBPT(tree,varargin)

% Removes all fields that have been introduced by any pruning function,
% returning a "clean" or "depruned" BPT
% Removes all fields that have been introduced by any function is specified
% by second string argument 'all'

if ~isprunedtree(tree) && nargin == 1
    
    disp('no ''pruning'' field in input tree. Output tree is the same as input tree.')
    
elseif nargin == 2
    
    opt = varargin{1};
    
    if ~strcmp(opt,'all')
        
        error('string option must be ''all'' to remove all additional fields.')
        
    end
    
    initialFields = {'label','descriptors','nodeinfo','construction'};
    allFieldsinTree = fieldnames(tree(1));
    
    if ~isequal(allFieldsinTree,fieldnames(tree(end)))
        error('Field names are not consistent along the tree structure array.')
    end
    
    strcomparison = false(size(allFieldsinTree));
    for i=1:length(strcomparison)
        strcomparison(i) = any(strcmp(allFieldsinTree{i},initialFields));
    end
    
    Fieldstoremove = allFieldsinTree(~strcomparison);
    tree = rmfield(tree,Fieldstoremove);
    
else
    
    error('Wrong number of input arguments.')
    
end