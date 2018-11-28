function BPT = getBPTenergy(BPT,data,initialseg,handleGofF,handleRegu)

% Tnew = getBPTenergy(Told,data,initialseg,GofF,Regu)
% This function adds a new field 'energy' to the BPT and retrieve for each
% node of the tree its two term energy, composed of the goodness-of-fit
% term 'GofF' and the regularization term 'Regu'. Both of them have to be
% handles, where GofF takes the pixel values of the region as an argument,
% and Regu takes a binary image representing the region.

fprintf('Retrieving energy for each node\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

clearvars -global
global tree
global im
global initsegmap

im = reshape(data,size(data,1)*size(data,2),size(data,3));
tree = BPT;
initsegmap = initialseg;

clearvars BPT data initialseg

Nbnodes = length(tree);
Nbleaves = (Nbnodes+1)/2;
Nbmergings = Nbnodes-Nbleaves;

for ii=1:Nbnodes
    tree(ii).energy = energystruct;
    
    % retrieve region
    Rleaves = tree(ii).nodeinfo.leaves;
    if isempty(Rleaves)
        Rleaves = ii;
    end
    R = ismember(initsegmap,Rleaves);
    valPixR = im(R(:),:);
    
    GofF = handleGofF(valPixR,ii); % get goodness-of-fit value
    Regu = handleRegu(R,ii); % get regularization term
    tree(ii).energy.properenergy = [GofF,Regu];
    
    % displaying some information to the user (function may take a while)
    switch ii
        case Nbleaves
            fprintf('Leaves: done\n')
        case Nbleaves + round(Nbmergings/4)
            fprintf('25%% of merged nodes: done\n')
        case Nbleaves + round(Nbmergings/2)
            fprintf('50%% of merged nodes: done\n')
        case Nbleaves + round(3*Nbmergings/4)
            fprintf('75%% of merged nodes: done\n')
    end
    
end

BPT = tree;
clearvars -global

function s = energystruct
s.properenergy = 0;
s.optimalenergy = 0;
s.optimalcut = 0;
s.childrenenergy = 0;