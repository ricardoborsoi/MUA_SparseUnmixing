function out = rescale(in,a)

[x y ~] = size(in);

m = min(min(in));
M = max(max(in));

m = repmat(m,[x,y,1]);
M = repmat(M,[x,y,1]);

out = a*(in-m)./(M-m);
if a == 1
    out = double(out);
end
if a == 255
    out = uint8(out);
end