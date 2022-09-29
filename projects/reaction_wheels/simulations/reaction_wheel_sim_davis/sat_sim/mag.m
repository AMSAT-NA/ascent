function y = mag( x, dim )%MAG Sum of squares.%   For vectors, MAG(X) is the sum of the squares of the %   elements of X. For matrices, MAG(X) is a row vector %   with the sum of the squares over each column. For N-D %   arrays, MAG(X) operates along the first non-singleton %   dimension.%%   MAG(X,DIM) operates along the dimension DIM. %%   Example: If X = [0 1 2%                    3 4 5]%%   then mag(X,1) is [3 17 29] and mag(X,2) is [5 50]'.%%   See also SUM and RMS.if nargin==1,	y = sum( x.*conj(x) );else	y = sum( x.*conj(x), dim );end