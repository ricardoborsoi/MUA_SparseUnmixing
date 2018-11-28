function segmapfalsecolor = displaysegmentationfalsecolors(segmap,imrgb)

[m,n,p] = size(imrgb);

if p == 1
    imrgb = repmat(imrgb,[1 1 3]);
    p = 3;
end

segmap = reshape(segmap,m*n,1);
segmapfalsecolor = zeros(m*n,p);

R0 = imrgb(:,:,1);
G0 = imrgb(:,:,2);
B0 = imrgb(:,:,3);

seglabels = unique(segmap(:));

for i=1:length(seglabels);
    ind = segmap==seglabels(i);
    R = R0(ind);
    G = G0(ind);
    B = B0(ind);
    r = mean(R);
    g = mean(G);
    b = mean(B);
    segmapfalsecolor(ind,:) = repmat([r g b],sum(ind),1);
end

segmapfalsecolor = reshape(segmapfalsecolor,m,n,p);    