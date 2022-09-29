function y = quat2axang(q)
%QUAT2AXANG Convert quaternion to axis-angle form

% see https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

y = [ q(2:4)/norm(q(2:4)) 2*acos(q(1)) ];