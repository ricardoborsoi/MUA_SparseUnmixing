function N = neigh4(A,i,j)

N(:,1) = A(i-1,j,:);
N(:,2) = A(i,j-1,:);
N(:,3) = A(i,j+1,:);
N(:,4) = A(i+1,j,:);