function [F,nbcall] = priority_size(merging_order,nbcall)

global tree
global I

if ~strcmp(merging_order,'stop')
    
    alpha = 0.15;
    
    Nregions = length(merging_order);
    mean_size_region = size(I,1)/Nregions;
    
    thresh = round(alpha*mean_size_region)  ;
    % the merging priority is given to all the regions whose size is below
    % thresh
    
end

persistent regsize
persistent regmax

if isempty(nbcall)
    nbcall = 1;
    regsize = zeros(size(tree));
    regmax = max(merging_order);
    s = [tree(1:regmax).descriptors];
    s = [s.size];
    regsize(1:regmax) = s;
elseif strcmp(merging_order,'stop')
    clear regsize
    clear regmax
    return
else
    nbcall = nbcall + 1;
    s = [tree(regmax:max(merging_order)).descriptors];
    s = [s.size];
    regsize(regmax:max(merging_order)) = s;
    regmax = max(merging_order);
end

F = regsize(merging_order)<thresh;
% F is a flag binary vector where 1 means that the concerned regions gets
% the merging priority