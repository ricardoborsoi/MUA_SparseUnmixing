function updatetree(newlabel)

global tree

child = tree(newlabel).nodeinfo.children;
lbli = child(1); lblj = child(2);
N = tree(newlabel).construction.neighbors;
D = tree(newlabel).construction.alldist;

for i=1:length(N)
    n = N(i);
    d = D(i);
    
    ind = ismember(tree(n).construction.neighbors,[lbli,lblj]);
    tree(n).construction.neighbors(ind) = [];
    tree(n).construction.alldist(ind) = [];
    
    tree(n).construction.neighbors(end+1) = newlabel;
    tree(n).construction.alldist(end+1) = d;
    
    if d<tree(n).construction.dist
        tree(n).nodeinfo.sibling = newlabel;
    end
    
    if ismember([lbli,lblj],tree(n).nodeinfo.sibling)
        ind = find(tree(n).construction.alldist==min(tree(n).construction.alldist),1,'first');
        tree(n).nodeinfo.sibling = tree(n).construction.alldist(ind);
    end
    
end