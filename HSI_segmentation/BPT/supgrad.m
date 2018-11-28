function [grad,map] = supgrad(in)

se = strel('square',3);
gr = zeros(size(in));
d = size(in,3);
for i=1:d
    gr(:,:,i) = imdilate(in(:,:,i),se)-imerode(in(:,:,i),se);
    gr(:,:,i) = rescale(gr(:,:,i),1);
end

grad = zeros(size(gr,1), size(gr,2));
for i=1:size(gr,1)
    for j=1:size(gr,2)
        grad(i,j) = max(gr(i,j,:));
    end
end

map = double(watershed(grad));