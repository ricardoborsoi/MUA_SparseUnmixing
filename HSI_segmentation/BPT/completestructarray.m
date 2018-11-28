function completestructarray

global tree

tree(end).nodeinfo.height = 1;
tree(end).nodeinfo.branch = [];

child = tree(end).nodeinfo.children;

while ~isempty(child);
    
    parentlabel = tree(child(1)).nodeinfo.parent;
    
    parents = cat(2,parentlabel,tree(parentlabel).nodeinfo.branch);
    tree(child(1)).nodeinfo.branch = parents;
    tree(child(1)).nodeinfo.height = length(parents)+1;
    tree(child(2)).nodeinfo.branch = parents;
    tree(child(2)).nodeinfo.height = length(parents)+1;
    
    if ~isempty(tree(child(1)).nodeinfo.children);
        child = cat(2,child,tree(child(1)).nodeinfo.children);
    end
    
    if ~isempty(tree(child(2)).nodeinfo.children);
        child = cat(2,child,tree(child(2)).nodeinfo.children);
    end
    
    child([1 2]) = [];
    
end

% NodeInfo = [tree.nodeinfo];
% Alllvl = [NodeInfo.height];
% Maxlvl = max(Alllvl);
% 
% for i=1:length(tree)
%     tree(i).nodeinfo.Nheight = tree(i).nodeinfo.height/Maxlvl;
% end