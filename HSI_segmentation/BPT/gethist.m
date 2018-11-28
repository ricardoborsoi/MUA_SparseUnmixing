function hist = gethist(pixls,Nbin)

S = size(pixls,2);
hist = histc(pixls,0:1/Nbin:1,2)/S;