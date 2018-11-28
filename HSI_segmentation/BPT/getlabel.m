function lbl = getlabel(map,Zw,flag)

m = size(map,1);
m0 = m-2;

yw = ceil(Zw/m0);
xw = Zw-m0*(yw-1);

if flag == 4
    N = neigh4(map,xw+1,yw+1);
else
    N = neigh8(map,xw+1,yw+1);
    if ~any(N)
        map(1,:) = [];
        map(end,:) = [];
        map(:,1) = [];
        map(:,end) = [];
        s = 2;
        while ~any(N)
            newmap = padarray(map,[s s],'symmetric');
            N = neigh8(newmap,xw+s,yw+s,s);
            s = s+1;
        end
    end
end
N(N==0) = [];
lbl = unique(N);