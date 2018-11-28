function N = neigh8(A,i,j,varargin)

if nargin == 3
    s = 1;
elseif nargin == 4
    s = varargin{1};
else
    error('wrong number of input argument')
end

N = A(i-s:i+s,j-s:j+s);
N((2*s+1)*s+s+1) = [];