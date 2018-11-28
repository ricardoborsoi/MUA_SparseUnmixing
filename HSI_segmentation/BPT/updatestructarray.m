function updatestructarray(handleO,handlemerging,varargin)

fprintf('Construction of the tree\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

global tree

D = [tree.construction];
D = [D.dist];
[D,merging_order] = sort(D,'ascend');
S = [tree.nodeinfo];
S = [S.sibling];
S = S(merging_order);
N_iteration = length(S);

if nargin==3
    handlepriority = varargin{1};
    nbcall = [];
end
F = false(size(merging_order));

% looping until there are two regions remaining
for i=1:N_iteration-2

    if i==round(N_iteration/4)
        fprintf('merging: 25%%\n')
    elseif i==round(N_iteration/2)
        fprintf('merging: 50%%\n')
    elseif i==round(3*N_iteration/4)
        fprintf('merging: 75%%\n')
    end
    
    % Step 1: retrieve the two regions which are to merge

    if nargin==2
        merging_id = 1;
    else
        if nargin==3&&all(F==false)
            [F,nbcall] = handlepriority(merging_order,nbcall);
        end
        merging_id = find(F,1,'first');
        if isempty(merging_id)
            merging_id = 1;
        end
    end
    
    node_i = tree(merging_order(merging_id));
    node_j = tree(S(merging_id));
    
    node_j.nodeinfo.sibling = node_i.label;
    node_j.construction.dist = node_i.construction.dist;
    tree(S(merging_id)) = node_j;
    % update information on sibling node (changes something only if
    % priority function is holding
   
    S(merging_id) = [];
    merging_order(merging_id) = [];
    D(merging_id) = [];
    F(merging_id) = [];

    % Step 2: Create the new structure in tree
    newlabel = N_iteration+i;
    tree(newlabel).label = newlabel;
    tree(newlabel).descriptors = handlemerging(node_i,node_j);
    updatenewnode(node_i.label,node_j.label,newlabel,handleO);
    tree(newlabel).nodeinfo.iteration = i;
    
    index = find(merging_order==node_j.label,1,'first');
    merging_order(index) = newlabel;
    S(index) = tree(newlabel).nodeinfo.sibling;
    D(index) = tree(newlabel).construction.dist;
    F(index) = false;
    
    % Step 3: update tree
    
    child = tree(newlabel).nodeinfo.children;
    lbli = child(1); lblj = child(2);
    Nei = tree(newlabel).construction.neighbors;
    Dist = tree(newlabel).construction.alldist;
    
    for j=1:length(Nei)
        n = Nei(j);
        d = Dist(j);
        
        ind = tree(n).construction.neighbors==lbli|tree(n).construction.neighbors==lblj;
        tree(n).construction.neighbors(ind) = [];
        tree(n).construction.alldist(ind) = [];
        
        tree(n).construction.neighbors(end+1) = newlabel;
        tree(n).construction.alldist(end+1) = d;
        
        if d<tree(n).construction.dist
            tree(n).construction.dist = d;
            tree(n).nodeinfo.sibling = newlabel;
            ind = merging_order==n;
            S(ind) = newlabel;
            D(ind) = d;
        end
        
        if any(tree(n).nodeinfo.sibling==[lbli,lblj])
            ind = find(min(tree(n).construction.alldist)==tree(n).construction.alldist,1,'first');
            tree(n).construction.dist = tree(n).construction.alldist(ind);
            tree(n).nodeinfo.sibling = tree(n).construction.neighbors(ind);
            ind = merging_order==n;
            S(ind) = tree(n).nodeinfo.sibling;
            D(ind) = tree(n).construction.dist;
        end
        
    end
    
    % Step 4: update merging order
    % If priority function, update only when no regions with priority
    % remain
    if all(F==false)
        [D,order] = sort(D,'ascend');
        merging_order = merging_order(order);
        S = S(order);
    end
    
end

% fusing the two last regions
node_i = tree(merging_order(1));
node_j = tree(S(1));
newlabel = length(tree);
tree(end).label = newlabel;
tree(end).descriptors = handlemerging(node_i,node_j);
updatenewnode(node_i.label,node_j.label,newlabel,handleO);
tree(end).nodeinfo.parent = [];
tree(end).nodeinfo.sibling = [];
tree(end).nodeinfo.issiblingleaf = false;
tree(end).nodeinfo.iteration = N_iteration-1;

fprintf('merging: 100%%\n')

% clearing persistent variables in priority function
if nargin==3
    handlepriority('stop',nbcall);
end