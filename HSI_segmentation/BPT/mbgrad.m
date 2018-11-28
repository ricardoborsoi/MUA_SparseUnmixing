function [grad,map] = mbgrad(in)

grad = zeros(size(in,1),size(in,2));
E = padarray(in, [1 1 0], 'symmetric');
[m n] = size(E(:,:,1));

for i=2:m-1
    for j=2:n-1
        N1 = E(i-1:i+1,j-1:j+1,:);
        N1 = reshape(N1,size(N1,1)*size(N1,2),size(N1,3));
        N1(5,:) = [];
        N1 = N1';
        x = squeeze(E(i,j,:));
        d1 = N1-repmat(x,1,size(N1,2));
        d2 = d1.^2;
        d3 = sum(d2);
        D = d3.^0.5;
        grad(i-1,j-1) = max(D)-min(D);
    end
end

map = double(watershed(grad));