function tree = getpersistentintervals(tree,varargin)

if nargin == 1
    display = false;
elseif nargin == 2
    opt = varargin{1};
    if strcmp(opt,'display')
        display = true;
    else
        error('Second input argument must be ''display''.')
    end
else
    error('Wrong number of input arguments.')
end

Nbnodes = length(tree);

% Get appearance value lambda+
for ii=1:Nbnodes

    properenergy = tree(ii).energy.properenergy;
    
    if isempty(tree(ii).nodeinfo.children)
        
        tree(ii).energy.lambdaplus = 0;
        tree(ii).energy.partialenergy = {[0,properenergy(2),properenergy(1)]};
        
    else
        
        nodeproperenergy = {[0,properenergy(2),properenergy(1)]};
        children = tree(ii).nodeinfo.children;
        child1partialenergy = tree(children(1)).energy.partialenergy;
        child2partialenergy = tree(children(2)).energy.partialenergy;
        childrenpartialenergy = addPAfunctions(child1partialenergy,child2partialenergy);
        [partialenergy,lambda_plus] = findpartialenergy(nodeproperenergy,childrenpartialenergy);
        tree(ii).energy.partialenergy = partialenergy;
        tree(ii).energy.lambdaplus= lambda_plus;
        
        if ii==Nbnodes && display == true
            drawfunctions(partialenergy,'root partial energy');
        end
        
    end
end

% get disappearance value lambda-
for jj=1:Nbnodes-1
    branch_jj = tree(jj).nodeinfo.branch;
    branchEnergy = [tree(branch_jj).energy];
    branch_lambdaplus = [branchEnergy.lambdaplus];
    tree(jj).energy.lambdamoins = min(branch_lambdaplus);
end
tree(end).energy.lambdamoins = Inf;

function f_sum = addPAfunctions(f_1,f_2)

% sum of two piecewise affine functions

Nb_breakpoints_1 = length(f_1);
Nb_breakpoints_2 = length(f_2);

breakpoints_1 = zeros(1,Nb_breakpoints_1); breakpoints_2 = zeros(1,Nb_breakpoints_2);
slope_1 = zeros(1,Nb_breakpoints_1); slope_2 = zeros(1,Nb_breakpoints_2);
offset_1 = zeros(1,Nb_breakpoints_1); offset_2 = zeros(1,Nb_breakpoints_2);

% Retrieve breakpoint, slope and offset of each segment of PA function #1
for ii=1:Nb_breakpoints_1
    breakpoints_1(ii) = f_1{ii}(1);
    slope_1(ii) = f_1{ii}(2);
    offset_1(ii) = f_1{ii}(3);
end

% Retrieve breakpoint, slope and offset of each segment of PA function #2
for jj=1:Nb_breakpoints_2
    breakpoints_2(jj) = f_2{jj}(1);
    slope_2(jj) = f_2{jj}(2);
    offset_2(jj) = f_2{jj}(3);
end

% sort all breakpoints in ascending order
all_breakpoints = cat(2,breakpoints_1,breakpoints_2);
sort_breakpoints = sort(all_breakpoints,'ascend');
sort_breakpoints = unique(sort_breakpoints);
Nb_sort_breakpoints = length(sort_breakpoints);

slope_sum = zeros(1,Nb_sort_breakpoints);
offset_sum = zeros(1,Nb_sort_breakpoints);
f_sum = cell(1,Nb_sort_breakpoints);

for kk=1:Nb_sort_breakpoints
    
    current_BP = sort_breakpoints(kk);
    
    idx_1 = find(breakpoints_1<=current_BP,1,'last');
    if isempty(idx_1)
        idx_1 = length(breakpoints_1);
    end
    current_slope_1 = slope_1(idx_1);
    current_offset_1 = offset_1(idx_1);
    
    idx_2 = find(breakpoints_2<=current_BP,1,'last');
    if isempty(idx_2)
        idx_2 = length(breakpoints_2);
    end
    current_slope_2 = slope_2(idx_2);
    current_offset_2 = offset_2(idx_2);
    
    slope_sum(kk) = current_slope_1+current_slope_2;
    offset_sum(kk) = current_offset_1+current_offset_2;
    
    f_sum{kk} = [sort_breakpoints(kk),slope_sum(kk),offset_sum(kk)];
    
end
[~,gaps] = findgaps(f_sum);
if any(gaps>1e-4)
    disp('fuck, another discontinuity...')
end

function [partialenergy,intersect_value] = findpartialenergy(f_1,f_2)

% Intersection between affine f_1 and piecewise affine f_2 functions

slope_node = f_1{1}(2);
offset_node = f_1{1}(3);

Nb_breakpoints = length(f_2);
breakpoints_children = zeros(1,Nb_breakpoints);
slope_children = zeros(1,Nb_breakpoints);
offset_children = zeros(1,Nb_breakpoints);

for ii=1:Nb_breakpoints
        breakpoints_children(ii) = f_2{ii}(1);
        slope_children(ii) = f_2{ii}(2);
        offset_children(ii) = f_2{ii}(3);
end

if slope_node == slope_children(end) && offset_node == offset_children(end)
    % see Laurent Guigues' phD, page 165
    intersect_value = breakpoints_children(end);
    
elseif slope_node == slope_children(end) && offset_node > offset_children(end)
    % see Laurent Guigues' phD, page 165
    intersect_value = +Inf;
    
else % Regularization term is stricly sub-additive
    for jj=Nb_breakpoints:-1:1
        
        intersect_val = (offset_node-offset_children(jj))/(slope_children(jj)-slope_node);
        if intersect_val < 0 || isnan(intersect_val)
            intersect_val = 0;
        end
        
        if jj==Nb_breakpoints
            if intersect_val>=breakpoints_children(jj)
                intersect_value = intersect_val;
                break;
            end
        elseif intersect_val<breakpoints_children(jj+1)&&intersect_val>=breakpoints_children(jj)
            intersect_value = intersect_val;
            break;
        end
        
    end
    
end

idx_children = find(breakpoints_children<=intersect_value,1,'last');
partialenergy = cell(1,idx_children+1);
for kk=1:idx_children
    partialenergy{kk} = [breakpoints_children(kk),slope_children(kk),offset_children(kk)];
end
if isinf(intersect_value)
    partialenergy(end) = [];
else
    partialenergy{end} = [intersect_value,slope_node,offset_node];
end

function [gaps_pos,gaps_val] = findgaps(f)

Nb_breakpoints = length(f);

breakpoints = zeros(1,Nb_breakpoints);
slope = zeros(1,Nb_breakpoints);
offset = zeros(1,Nb_breakpoints);
gaps_pos = zeros(1,Nb_breakpoints);
gaps_val = zeros(1,Nb_breakpoints);

for ii=1:Nb_breakpoints
    breakpoints(ii) = f{ii}(1);
    slope(ii) = f{ii}(2);
    offset(ii) = f{ii}(3);
end

% check for gaps at breakpoints
for ll=2:Nb_breakpoints
    left_value = f{ll-1}(3) + f{ll}(1)*f{ll-1}(2);
    right_value = f{ll}(3) + f{ll}(1)*f{ll}(2);
    diff_value = right_value-left_value;
    gaps_val(ll) = diff_value;
    gaps_pos(ll) = abs(diff_value/left_value);
%    gaps_pos(ll) = ~isequal(diff_value,0);
end

function drawfunctions(varargin)

drawfunct = figure;
colors = {'b','r','k','c','m'};


cmpt = 1;
for ii=1:2:length(varargin)
    f_all{cmpt,:} = varargin{ii};
    legend_i{cmpt,:} = varargin{ii+1};
    cmpt = cmpt+1;
end

for jj=1:length(f_all)
    
    % extract one PA function from the cell array
    f_i = f_all{jj};
    % Initialize breakpoints, slope and offset arrays
    Nb_breakpoints_i = length(f_i);
    breakpoints_i = zeros(1,Nb_breakpoints_i);
    slope_i = zeros(1,Nb_breakpoints_i);
    offset_i = zeros(1,Nb_breakpoints_i);
    
    % Fill breakpoint, slope and offset arrays
    for kk=1:Nb_breakpoints_i
        breakpoints_i(kk) = f_i{kk}(1);
        slope_i(kk) = f_i{kk}(2);
        offset_i(kk) = f_i{kk}(3);
    end
    
    x_i = zeros(20,Nb_breakpoints_i);
    y_i = zeros(size(x_i));
    breakpoints_pos = cat(2,0,repmat(breakpoints_i(2:end),[1,2]));
    breakpoints_pos = sort(breakpoints_pos,'ascend');
    breakpoints_val = zeros(size(breakpoints_pos));
    
    % Get x-axis and y-axis values
    for ll=1:Nb_breakpoints_i-1
        x_i(:,ll) = linspace(breakpoints_i(ll),breakpoints_i(ll+1),20);
        y_i(:,ll) = offset_i(ll)+x_i(:,ll)*slope_i(ll);
        breakpoints_val(2*ll-1) = offset_i(ll)+breakpoints_i(ll)*slope_i(ll);
        breakpoints_val(2*ll) = offset_i(ll)+breakpoints_i(ll+1)*slope_i(ll);
    end
    x_i(:,end) = linspace(breakpoints_i(end),3*breakpoints_i(end),20);
    y_i(:,end) = offset_i(end)+x_i(:,end)*slope_i(end);
    breakpoints_val(end) = offset_i(end)+breakpoints_i(end)*slope_i(end);
    
    figure(drawfunct)
    plot(x_i(:),y_i(:),'Color',colors{jj});
    hold on
    plot(breakpoints_pos,breakpoints_val,'Color',colors{jj},...
        'Linestyle','none','Marker','x','MarkerSize',10);
    
end

newlegend_i = cell(2*size(legend_i,1),1);
for mm=1:length(legend_i)
    newlegend_i{2*mm-1} = legend_i{mm};
    newlegend_i{2*mm} = ['breakpoints ',legend_i{mm}];
end

legend(newlegend_i)