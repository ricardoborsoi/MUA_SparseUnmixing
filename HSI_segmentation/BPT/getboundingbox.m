function B = getboundingbox(box,S)

% Gets bounding box for a region

bb = box.BoundingBox;

xmin = max(1,floor(bb(2)));
ymin = max(1,floor(bb(1)));
xmax = min(S(1),ceil(bb(2))+bb(4));
ymax = min(S(2),ceil(bb(1))+bb(3));

B = [xmin xmax ymin ymax];
