function nweight = normweight(weight)

% Normalizes the weight matrix so each band has its weights summed to 1

[m n p] = size(weight);
W = squeeze(sum(sum(weight)));
W = reshape(W,1,1,p);
sumweight = repmat(W,[m n]);
nweight = weight./sumweight;