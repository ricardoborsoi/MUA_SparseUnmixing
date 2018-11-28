function Regu = Regu_boundary(R,~)

Bound = bwmorph(R,'remove'); 
Perim = sum(Bound(:));
Regu = Perim/2;
