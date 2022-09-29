clear
close all

%SIM_MAIN Simulate a satellite in inertial space
%   C struct -- contains most sim initial parameters.
%   S struct -- sim variables passed between the sim and the controller.
%   D struct -- saves sim results for plots.
%   sim_plot.m plots sim results using the C and D structs.
%   sim_3d_plot.m draws the satellite orientations as the sim runs.

t0 = 0;										% start controller @ t0
c.sat = '3u';
c.cntrl = 'sm'; 							% control algorithm
c.failed_motor = [];

% Body constants
%
V = [ 1  1 -1 -1  1  1 -1 -1;				% vertices of unit cube
	  1 -1 -1  1  1 -1 -1  1;
	  1  1  1  1 -1 -1 -1 -1];

switch c.sat
	case 'kok'								% Kok satellite
		c.Jb = diag([4 4 3]);				% inertia matrix (Kg*m^2)
		V = 3*V.*[4 4 3]';					% dimensions
		c.Jw = 5e-4 * eye(4);				% wheel inertias (Kg*m^2)
		c.K = 12.34e-3;						% motor torque constant
		endT = 100;
		t1 = 50; 							% reference trajectory time (sec)
	
	case '3u'								% 3U satellite
		c.Jb = diag([3.33e-02 3.33e-02 6.67e-03]);% inertia matrix (Kg*m^2)
		V = 5*V.*[1 1 3]';					% dimensions
		c.Jw = 6.8682e-05 * eye(4);			% wheel inertias (Kg*m^2)
		c.K = 17e-3;						% motor torque constant
		endT = 10;
		t1 = 10; 							% reference trajectory time (sec)
	
	case '6u'								% 3U satellite
		c.Jb = diag([8.67e-02 6.67e-02 3.33e-02]);% inertia matrix (Kg*m^2)
		V = 5*V.*[1 2 3]';					% dimensions
		c.Jw = 6.8682e-05 * eye(4);			% wheel inertias (Kg*m^2)
		c.K = 17e-3;						% motor torque constant
		endT = 20;
		t1 = 20; 							% reference trajectory time (sec)
end


% direction cosines of wheel axes
%
c.L = [ 1 -1 -1  1;							% axis or wheels
	    1  1 -1 -1;
	    1 -1  1 -1 ] / sqrt(3);
c.L(:,c.failed_motor)=0;


% Body initial & final conditions
%
eul0 = [0 0 0] * pi/180;					% initial [yaw pitch roll] of body
eul1 = [30 20 40] * pi/180;					% final [yaw pitch roll] of body

rpm_b = -[0 0 0]'*0.1;    					% body RPM
s.omega_b = rpm_b*pi/30;					% body rate (rad/sec)

% Wheel initial conditions
%
rpm_w = 000;									% wheel RPM
s.omega_w = (pinv(c.L)*[1 0 0]')*pi/30*rpm_w;	% wheel rates (rad/sec)

c.q0 = eul2quat(eul0);
c.q1 = eul2quat(eul1);
s.q = quatdivide(c.q0,c.q1);					% q = q0/q1	 q1*q = q0

c.delT = endT/1000;
N = endT/c.delT;
t = (0:N)' * c.delT; c.t = t;


% disturbance torque on body
%
c.Td = 0.0002 * randn(3,1) .* ones(3,N+1);
c.Td(:, t<t(end)/2 | t>t(end)*2/3 ) = 0;


% Reference Trajectory S-curves
%
c.ref = ones(N+1,1);
tt = t(t<t1)/t1;
c.ref(t>=t0 & t<t0+t1) = (1 - cos( pi*tt ) ) / 2;
% c.ref(t>=t0 & t<t0+t1) = tt.^2;
% c.ref =  1 ./ ( 1+exp(-12*(t-(t0+t1/2))/t1) );
%	figure; plot( t, c.ref, 'b' ); grid; %return
	
	
d.q = zeros(N+1,4);								% q
d.omega_b = zeros(3,N+1);						% body rate (rad/sec)
d.omega_w = zeros(4,N+1);						% wheel rate (rad/sec)
d.Tw = zeros(3,N+1);							% wheel torque (Nm)

R = quat2dcm(c.q0)';							% R = quat2dcm(q0)'
sim_3d_plot( 'init', R, V, 0 );
controller( ['init_' c.cntrl], c, s );

Tw = [0 0 0]';
for i = 1:N+1
	
	if t(i)>=t0								% start control @ t=2 sec
		Tw = controller( ['run_' c.cntrl], c, s, c.ref(i) );
%		Tw = Tw .* ( 1+.5*sin( 2*pi*t(i)/2 ) );
%		Tw = Tw .* ( 1+.1*randn(3,1) );
	end
	
	% save data for plots after sim completes
	%
	d.q(i,:) = s.q;								% position error
	d.omega_b(:,i) = s.omega_b;					% body rate (rad/sec)	
	d.omega_w(:,i) = s.omega_w;					% wheel rates (rad/sec)
	d.Tw(:,i) = Tw;								% motor torque (Nm)
	
	
	% Calculate angular acceleration on body in body coordinates
	%
	hw = c.L*(c.Jw * s.omega_w);				% wheel angular momentum
	hb = c.Jb * s.omega_b;						% body angular momentum
 	Tb = Tw + c.Td(:,i) - cross(s.omega_b,hb+hw);		% body torque 
	omega_b_dot = c.Jb \ Tb;					% body acceleration
	s.omega_b = s.omega_b + omega_b_dot*c.delT;
	

	% integrate in body coordinates
	%
	omega_w_dot = pinv(c.L*c.Jw)*Tw;
	s.omega_w = s.omega_w - omega_w_dot*c.delT;
		
	s.q = s.q + qdot(s.q,s.omega_b)*c.delT;
	s.q = s.q / norm(s.q);						% reduce error of q
	
	
	% update 3-D plot every n steps. sim runs slower with small n.
	% 
	n = 2^5;
	if mod(i,n) == 0
		R = quat2dcm(quatmultiply(c.q1,s.q))';
		sim_3d_plot( 'update', R, V, 1e2*s.omega_b )
	end
end

R = quat2dcm(quatmultiply(c.q1,s.q))';
sim_3d_plot( 'update', R, V, 1e2*s.omega_b )		% final update

axang = quat2axang(s.q);
fprintf( 'Final position err @ t(end) = %.2f m deg\n', axang(4)*30e3/pi )


sim_plots( c, d )


function Tw = controller( op, c, s, ref )
% CONTROLLER -- drive s.q(2:4) and s.omega_b to zero with torque feedback.
%	c = controller( op, c, s ) where op = 'init_lqr' or 'init_sm'
%  Tw = controller( op, c, s ) where op = 'run_lqr' or 'run_sm'

	persistent K G A x


	switch op
		
		case 'init_lqr'					% init LQR Controller
			hw = c.L*(c.Jw * s.omega_w);
			A = [zeros(3) eye(3)/2; zeros(3) c.Jb\skew(hw)];
			B = [zeros(3); pinv(c.Jb)];
			C = eye(6);
			Q = eye(6);
			R = 100*eye(3);
			if strcmp( c.sat, 'kok' )
				K = lqr( A, B, Q, R/32 );
			else							% 3U & 6U bodies
				K = lqr( A, B, Q, R );
				K(:,1:3) = 4*K(:,1:3);
			end
			
		case 'run_lqr'					% LQR Controller
			Tw = -K * [ref*s.q(2:4)'; s.omega_b];


		case 'init_sm'					% init Sliding Mode Controller
			if strcmp( c.sat, 'kok' )
				K = 0.3;
				G = 0.1;
			else							% 3U body
				K = 0.7*10;
				G = 1;
			end
			
		case 'run_sm'					% Sliding Mode Controller
			hw = c.L*(c.Jw * s.omega_w);	% angular momentum of wheels
			hb = c.Jb * s.omega_b;			% angular momentum of body
			q_dot = qdot(s.q,s.omega_b)';
			Tw = cross(s.omega_b,hb+hw) - c.Jb*( K*q_dot(2:4) ...
				+ G*sgn(s.omega_b + K*s.q(2:4)'*ref) );

		case 'init_um'
			axang = quat2axang(s.q);
			angAx = axang(1:3)'*axang(4);
			if strcmp( c.sat, 'kok' )
				load umich_kok
			else							% 3U body
				load umich_3u
			end
			% 			Obs = [c.C; c.C*c.A];
% 			y0 = [ ( c.L \ angAx )'; zeros(5,4) ];
% 			x = Obs \ y0;
 			x = A \ [zeros(3,4); ( pinv(c.L)*angAx )' ];
			
		case 'run_um'
			hw = c.L*(c.Jw * s.omega_w);				% angular momentum of wheels
			hb = c.Jb * s.omega_b;						% angular momentum of body
			
			u = [[cross(s.omega_b,hb+hw)' 0];zeros(3,4)];
			x_dot = A*x;
			x = x + x_dot * c.delT;
		%	x(1,:) = ( c.L \ angAx )';
			
			Tw = c.L * x(3,:)' .* diag(c.Jb)/c.Jb(1) * c.K  ...
				+ cross(s.omega_b,hb+hw);	
	
	end
end


function y = sgn(x)
%SGN Sliding mode surface thickness function

gain = 4;				% high gain => faster settling; higher current
%y = sign(gain*x);					% option A (chatters)
%y = max( -1, min(1,gain*x) );		% option B (discontinuous torques)
%y = tanh(gain*x);					% option C (smooth)
y = gain*x; 						% option D (smooth)
end

