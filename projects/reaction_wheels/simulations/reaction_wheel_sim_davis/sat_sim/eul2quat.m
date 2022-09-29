function q = eul2quat(e)
%EUL2QUAT Convert Euler (yaw,pitch,roll) angles to quaternion

% see https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

c = cos(e/2);
s = sin(e/2);

q = [c(1).*c(2).*c(3) + s(1).*s(2).*s(3), ...
	 c(1).*c(2).*s(3) - s(1).*s(2).*c(3), ...
	 c(1).*s(2).*c(3) + s(1).*c(2).*s(3), ...
	 s(1).*c(2).*c(3) - c(1).*s(2).*s(3)];
