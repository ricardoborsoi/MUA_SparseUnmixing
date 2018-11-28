function tree = pruneBPTnbregions(tree,NbReg)

N = length(tree);
Nbleaves = (N+1)/2;
for i=1:N
    tree(i).pruning = 0;
end

NbReg = round(NbReg);

if NbReg>Nbleaves
    
    for i=1:Nbleaves
        tree(i).pruning = 1;
    end
    
elseif NbReg>0
    
    IterMax = tree(end).nodeinfo.iteration;
    tree = pruneBPTnbiterations(tree,IterMax-NbReg+1);
    
else
    error('desired number of regions must be a positive integer')    
end