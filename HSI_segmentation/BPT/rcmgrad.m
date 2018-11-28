function [grad,map] = rcmgrad(in)

grad = zeros(size(in,1),size(in,2));
E = padarray(in, [1 1 0], 'symmetric');
[m n] = size(E(:,:,1));

for i=2:m-1
    for j=2:n-1
        N2 = E(i-1:i+1,j-1:j+1,:);
        N2 = reshape(N2,size(N2,3),size(N2,1)^2);
        di = pdist(N2','seuclidean');
        di = unique(di);
        di = sort(di,'descend');
        if length(di)>=2
            val = di(2);
        else
            val = di;
        end
        grad(i-1,j-1) = val;
    end
end

map = double(watershed(grad));