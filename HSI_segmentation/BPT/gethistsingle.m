function hist = gethistsingle(lbl,im,Nbin)

global initsegmap

[m n d] = size(im);

M = max(max(im));
mi = min(min(im));
M = repmat(M,[m,n]);
mi = repmat(mi,[m,n]);
im = (im-mi)./(M-mi);

S = 3; % Search window 7x7
N = 1; % Neighborhood window 3x3

im2 = padarray(im,[S+N S+N 0],'symmetric');
[x y] = find(initsegmap == lbl);
N1 = im2(x+S:x+S+2*N,y+S:y+S+2*N,:);
% d = size(N1,3);
we = zeros(2*S+1,2*S+1,d);

for i=x+N:x+N+2*S
    for j=y+N:y+N+2*S
        if (i~=x+N+S)||(j~=y+N+S)
            N2 = im2(i-N:i+N,j-N:j+N,:);
            we(i-(x+N)+1,j-(y+N)+1,:) = weight(N1,N2);
        else
            we(i-(x+N)+1,j-(y+N)+1,:) = zeros(1,d);
        end
    end
end

we(S+1,S+1,:) = max3d(we);
nwe = normweight(we);
hist = buildhistsingle(nwe,im2(x+N:x+N+2*S,y+N:y+N+2*S,:),Nbin);