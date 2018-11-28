function tree = pruneBPTnbiterations(tree,NbIter)

N = length(tree);
Nbleaves = (N+1)/2;
for i=1:N
    tree(i).pruning = 0;
end

IterMax = tree(end).nodeinfo.iteration;
NbIter = round(NbIter);

if NbIter>=IterMax
    tree(end).pruning = 1;
else
    IterCut = NbIter;
    SetLeaves = 1:Nbleaves;
    
    SetNode = [];
    while ~isempty(SetLeaves)
        leaf = SetLeaves(1);
        p = tree(leaf).nodeinfo.branch;
        IterBranch = [tree(p).nodeinfo];
        IterBranch = [IterBranch.iteration];
        node = p(find(IterBranch<=IterCut,1,'last'));
        if isempty(node)
            node = leaf;
            SetLeaves(1) = [];
        else
            leaves = tree(node).nodeinfo.leaves;
            SetLeaves = setdiff(SetLeaves,leaves);
        end
        SetNode = cat(2,SetNode,node);
    end
    
    for i=1:length(SetNode)
        tree(SetNode(i)).pruning = 1;
    end
    
end

