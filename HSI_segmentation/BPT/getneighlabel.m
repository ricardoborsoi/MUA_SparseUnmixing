function neighlbl = getneighlabel(map,lbl)

bwmap = (map==lbl);
se = strel('arbitrary',[0 1 0;1 1 1;0 1 0]);
bwmap2 = imdilate(bwmap,se)-bwmap;
neighlbl = unique(double(map.*bwmap2))';
neighlbl(neighlbl==0) = [];