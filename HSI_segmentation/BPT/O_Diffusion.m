function dist = O_Diffusion(Hi,Hj)

% computes the diffusion distance bewteen the two sets of histograms Hi and
% Hj

L = 3;
d(:,:,1) = Hi-Hj;
m = size(d,2)-1;
sig = 1;
t = -3:6/m:3;
gker = 1/(sig*sqrt(2*pi))*exp(-t.^2/(2*sig^2));

v=[1:2:2*m+1];
for i=2:L+1
    C = conv2(1,gker,d(:,:,i-1));
    d(:,:,i) = C(:,v);
end

dist = sum(sum(sum(abs(d))));