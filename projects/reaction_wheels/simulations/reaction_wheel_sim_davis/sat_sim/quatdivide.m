function q = quatdivide( s, r )
%QUATDIVIDE Divide quaternions: Q = S/R.

% see https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

r(2:4) = -r(2:4);
q = quatmultiply( r, s );