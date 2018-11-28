function mmx = max3d(M)

% M is a 3-D matrix. Computes the maximum value for each band of M.

mmx = max(max(M,[],2),[],1);
