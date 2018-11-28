function tree = initstructarray(handleR,handleO)

fprintf('Initialization of the tree\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

global initsegmap

[m, n] = size(initsegmap);
N = max(initsegmap(:));
tree(2*N-1) = struct;
boxes = regionprops(initsegmap,'Boundingbox');

for i=1:N
    tree(i).label = i;
    tree(i).descriptors = handleR(i);
    bb = getboundingbox(boxes(i),[m n]);
    tree(i).descriptors.boundingbox = bb;
    tree(i).nodeinfo = initnodeinfo;
    tree(i).construction.neighbors = getneighlabel(initsegmap(bb(1):bb(2),bb(3):bb(4)),i);
    tree(i).construction.alldist = zeros(size(tree(i).construction.neighbors));
end
fprintf('Creating leaf structures: done\n')

for i=1:N
    neigh = tree(i).construction.neighbors;
    Nei = neigh(neigh>tree(i).label);
    D = zeros(size(Nei));
    for j=1:length(Nei)
        D(j) = handleO(tree(i).descriptors.model,tree(Nei(j)).descriptors.model);
        ind = tree(Nei(j)).construction.neighbors == i;
        tree(Nei(j)).construction.alldist(ind) = D(j);
    end
    tree(i).construction.alldist(neigh>tree(i).label) = D;
    tree(i).construction.dist = min(tree(i).construction.alldist);
    ind = find(tree(i).construction.alldist==min(tree(i).construction.alldist),1,'first');
    tree(i).nodeinfo.sibling = neigh(ind);
    tree(i).nodeinfo.leaves = [];
    tree(i).nodeinfo.iteration = 0;
end
fprintf('Updating leaf structures: done\n')