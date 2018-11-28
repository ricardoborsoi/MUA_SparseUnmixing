function T = buildBPT(inputdata)

% define global variables
clearvars -global
global initsegmap
global I
global tree

% extract input parameters from input structure
initsegmap = inputdata.initsegmap;
I = inputdata.data; I = reshape(I,size(I,1)*size(I,2),size(I,3));
R = inputdata.regionmodel;
O = inputdata.mergingcriterion;
spec_merging = inputdata.specmerging;
priority_size = inputdata.prioritysize;

% BPT initialization step
tic
tree = initstructarray(R,O);
toc

% BPT update step
tic
updatestructarray(O,spec_merging,priority_size);
toc

% BPT completion step
completestructarray

% output
T = tree;
clearvars -global