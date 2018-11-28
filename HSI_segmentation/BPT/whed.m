function mapout = whed(mapin,I)

med = getmedian(mapin,I);
mapout = mapin;
Xw = find(mapin==0);
mapin = padarray(mapin,[1 1],'symmetric');

for i=1:length(Xw)
    x = I(Xw(i),:);
    lbl = getlabel(mapin,Xw(i),8);
    distance = zeros(size(lbl));
    for j=1:length(lbl)
        y = med(lbl(j),:);
        distance(j) = (sum(x-y).^2).^0.5;
    end
    indice = find(distance == min(distance));
    mapout(Xw(i)) = lbl(indice(1));
end