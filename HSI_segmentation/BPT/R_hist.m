function D = R_hist(lbl)

% The region model of the region defined by pixels is its set of normalized
% histograms (one histogram per band), meaning that for each band, the 
% histogram bins sum up to 1.

global initsegmap
global I
pixels = getpixels(I,initsegmap,lbl)';

D = struct;
D.size = size(pixels,2);

M = max(I,[],1)';
M = repmat(M,[1,size(pixels,2)]);
m = min(I,[],1)';
m = repmat(m,[1,size(pixels,2)]);

Nbin = 64; %number of bins in the histograms (has to be an even number)

if size(pixels,2)>1
    pixels = ((pixels-m)./(M-m));
    D.model = gethist(pixels,Nbin);
else
    [row col] = size(initsegmap);
    dim = size(I,2);
    D.model = gethistsingle(lbl,reshape(I,row,col,dim),Nbin);
end