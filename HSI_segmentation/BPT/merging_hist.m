function D = merging_hist(node_i,node_j)

% Specifies the region model of the merging of two regions Ri and Rj. The
% region model is the set of histograms of the region.

Di = node_i.descriptors;
Dj = node_j.descriptors;

% field size
Si = Di.size;
Sj = Dj.size;
D.size = Si+Sj;

%field model
Ri = Di.model;
Rj = Dj.model;
D.model = (Si*Ri+Sj*Rj)/(Si+Sj);

%field boundingbox
bbi = Di.boundingbox;
bbj = Dj.boundingbox;
D.boundingbox = [min(bbi(1),bbj(1)) max(bbi(2),bbj(2)) ...
    min(bbi(3),bbj(3)) max(bbi(4),bbj(4))];