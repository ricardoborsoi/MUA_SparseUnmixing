function H = buildhistsingle(W,N,Nbin)

% Derives the histograms from the weight matrix W and the corresponding
% neighbordhood values N

[m n p] = size(N);
H = zeros(p,Nbin+1);
W = reshape(W,m*n,p)';
N = reshape(N,m*n,p)';
[~,bin] = histc(N,0:1/Nbin:1,2);

for i=1:p
    H(i,:) = accumarray(bin(i,:)',W(i,:)',[Nbin+1,1])';
end