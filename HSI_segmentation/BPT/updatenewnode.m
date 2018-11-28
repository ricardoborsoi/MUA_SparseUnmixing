function updatenewnode(lbli,lblj,newlabel,handleO)

global tree

% filling in the nodeinfo structure for the new node, and updating it for
% the children

tree(lbli).nodeinfo.parent = newlabel;
tree(lblj).nodeinfo.parent = newlabel;

tree(newlabel).nodeinfo.children = [lbli lblj];

li = tree(lbli).nodeinfo.leaves;
lj = tree(lblj).nodeinfo.leaves;
if isempty(li)
    li = lbli;
    tree(lblj).nodeinfo.issiblingleaf = true;
end
if isempty(lj)
    lj = lblj;
    tree(lbli).nodeinfo.issiblingleaf = true;
end
tree(newlabel).nodeinfo.leaves = [li lj];
tree(newlabel).nodeinfo.nbleaves = length([li lj]);
tree(newlabel).nodeinfo.issiblingleaf = false;

% filling in the construction structure for new node (except if root)

if newlabel~=length(tree)
    % field newneighbors
    neighbors_i = tree(lbli).construction.neighbors;
    neighbors_j = tree(lblj).construction.neighbors;
    newneighbors = sort(union(neighbors_i,neighbors_j));
    newneighbors(newneighbors==lbli|newneighbors==lblj) = [];
    tree(newlabel).construction.neighbors = newneighbors;
    
    % fields alldist
    tree(newlabel).construction.alldist = zeros(size(newneighbors));
    Rnewlabel = tree(newlabel).descriptors.model;
    
    for i=1:length(newneighbors)
        tree(newlabel).construction.alldist(i) = handleO(Rnewlabel,...
            tree(newneighbors(i)).descriptors.model);
    end
    tree(newlabel).construction.dist = min(tree(newlabel).construction.alldist);
    
    % field dist + field sibling in structure nodeinfo
    ind = find(tree(newlabel).construction.dist==tree(newlabel).construction.alldist,1,'first');
    tree(newlabel).nodeinfo.sibling = newneighbors(ind);
end
