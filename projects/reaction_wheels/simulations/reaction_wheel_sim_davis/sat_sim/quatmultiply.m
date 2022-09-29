function q = quatmultiply( r, s )
%QUATMULTIPLY The product of two quaternions: Q = R*S.

% see https://mathworld.wolfram.com/Quaternion.html

q = [ (r(:,1)*s(:,1) - r(:,2)*s(:,2) - r(:,3)*s(:,3) - r(:,4)*s(:,4)), ...
	  (r(:,1)*s(:,2) + r(:,2)*s(:,1) + r(:,3)*s(:,4) - r(:,4)*s(:,3)), ...
	  (r(:,1)*s(:,3) - r(:,2)*s(:,4) + r(:,3)*s(:,1) + r(:,4)*s(:,2)), ...
	  (r(:,1)*s(:,4) + r(:,2)*s(:,3) - r(:,3)*s(:,2) + r(:,4)*s(:,1)) ];


