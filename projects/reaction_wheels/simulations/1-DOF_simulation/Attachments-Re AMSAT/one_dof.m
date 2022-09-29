% Cube dynamics
%	Really simple 1-DOF for now

clear
close all


% unit vectors of wheel axes (pyramid configuration)
%
V = [ 1  1 -1 -1;
	  1 -1  1 -1;
	  1 -1 -1  1 ] / sqrt(3);


% Momentum relations between body and wheels 
%
A = V;								%  body = A * wheel
B = pinv(A);						% wheel = B * body
% A*B = eye(3)						% check 


% Use largest inertia of satallite, Ixx (Kg*m^2)
%
Ib = 6.944e-05; 					% 0U5
Ib = 2.222e-03; 					% 1U
Ib = 5.417e-03; 					% 1U5
Ib = 3.333e-02; 					% 3U
Ib = 5.333e-01; 					% 12U

Tm = 10e-3;							% rated torque of motor (Nm)
Iw = 73.9e-6;						% wheel inertia (Kg*m^2)
									%  = 739 (g*cm^2)

% run motors at max |torque|
%
temp = B*[1 0 0]';
Tw = temp / max(abs(temp)) * Tm;	% max |torque| to any wheel = Tm

Tb = A * Tw;						% torque on body (Nm) from wheels  


% time to spin up 1 wheel from 0 RPM to 6000 RPM
%
rpm = 6e3;							% (rev/min)
Vw = pi * rpm/30;					% wheel velocity (rad/sec)
delT = Iw * Vw / Tm;				% (sec)


% calculate body velocity after delT secs 
% 
Vb = Tb(1) * delT / Ib;				% body velocity (rad/sec)


% Check that momentum of wheels & body are equal
%
%		Mw is momentum vector from all four wheels at 6000 RPM
%
Mw  = A * Vw * Iw*[1 1 -1 -1]';		% total wheel momentum (Nms)
Mb = Ib * Vb;						% body momentum (Nms)
% Mb - Mw(1)

fprintf( 'Time = %5.3f sec   Vb = %4.1f °/s   RPM = %4.1f\n', ...
	delT, 180/pi*Vb, 30/pi*Vb  )

