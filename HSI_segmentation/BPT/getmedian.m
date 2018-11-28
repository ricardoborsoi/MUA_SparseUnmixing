function med = getmedian(map,I)

med = zeros(max(map(:)),size(I,2));
for i=1:max(map(:))
    pix = getpixels(I,map,i);
    med(i,:) = median(pix);
end