function D = R_mean(lbl)

% Creates a structure. The region model of the set of pixels pixls is its
% mean spectrum

global initsegmap
global I
pixels = getpixels(I,initsegmap,lbl)';

D = struct;
D.size = size(pixels,2);
D.model = mean(pixels,2);