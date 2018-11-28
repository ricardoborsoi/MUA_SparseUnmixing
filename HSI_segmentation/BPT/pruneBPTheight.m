function tree = pruneBPTheight(tree,h)

N = length(tree);
Nbleaves = (N+1)/2;
for i=1:N
    tree(i).pruning = 0;
end

H = [tree.nodeinfo];
H = [H.height];

ind = false(1,N);
ind(1:Nbleaves) = H(1:Nbleaves)<=h;
ind(Nbleaves+1:end) = H(Nbleaves+1:end)==h;

node = find(ind);
for i=1:length(node)
    tree(node(i)).pruning = 1;
end