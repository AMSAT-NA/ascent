function eul = quat2eul( quat )
%QUAT2EUL Convert quaternions to Euler angles

% see https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

q1 = quat(:,1);
q2 = quat(:,2);
q3 = quat(:,3);
q4 = quat(:,4);

eul = zeros(size(quat,1), 3 );

sinel = -2*(q2.*q4-q1.*q3);
sinel(sinel >  1) =  1;
sinel(sinel < -1) = -1;


eul = [ atan2( 2*(q2.*q3+q1.*q4), q1.^2 + q2.^2 - q3.^2 - q4.^2 ), ...
		asin( sinel ), ...
		atan2( 2*(q3.*q4+q1.*q2), q1.^2 - q2.^2 - q3.^2 + q4.^2 ) ];
