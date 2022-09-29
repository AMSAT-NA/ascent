function sim_plots( c, d )
%SIM_PLOTS Miscellanious simulation plots

fig_info = ['sat = ' c.sat, '  cntrl = ', c.cntrl, ...
	'   failed motor = ',  num2str(c.failed_motor)];

figure( 'name', fig_info )
	plot( c.t, 2*acos(d.q(:,1)), 'b',						...
		  c.t, sqrt(mag(d.omega_b))*30/pi, 'r',	...
		  c.t, sqrt(mag(d.omega_w))*30e-3/pi, 'g'); hold on
	plot( c.t, sqrt(mag(d.Tw))/c.K, 'k' )
	legend('Position err (rad)', 'Body (rpm)', 'Wheel (Krpm)', 'Current (amps)' )
	grid
	zoomrb

	
eul = quat2eul(quatmultiply(c.q1,d.q))*180/pi;
figure( 'name', fig_info, 'NumberTitle', 'off' )
	annotation('textbox', [0, .53, .5, 0], 'string', fig_info, 'LineStyle', 'none' );
	subplot(2,1,1)
	plot( c.t, eul(:,1), 'r', c.t, eul(:,2), 'g', c.t, eul(:,3), 'b' )
	ylabel( 'Angle (deg)' ); xlabel( 'Time (sec)' );
	title( 'Satellite Euler angles' )
	legend('Yaw', 'Pitch', 'Roll' )
	grid
	
	subplot(2,1,2)
	rpm = d.omega_b*1e3;
	plot( c.t, rpm(1,:), 'r', c.t, rpm(2,:), 'g', c.t, rpm(3,:), 'b' )
	ylabel( 'Body rate (mrad/sec)' ); xlabel( 'Time (sec)' );
	title( 'Satellite Angular Velocity' )
	legend('Yaw', 'Pitch', 'Roll' )
	grid
	zoomrb
	
	
figure( 'name', fig_info, 'NumberTitle', 'off' )
	annotation('textbox', [0, .53, .5, 0], 'string', fig_info, 'LineStyle', 'none' );
	subplot(2,1,1)
	tw = pinv(c.L)*d.Tw*1e3;
	plot( c.t, tw(1,:), 'r', c.t, tw(2,:), 'g', c.t, tw(3,:), 'b', c.t, tw(4,:), 'k' )
	ylabel( 'Torque (mNm)' ); xlabel( 'Time (sec)' );
	title( 'Reaction Wheel Torgues' )
	legend('Wheel 1', 'Wheel 2', 'Wheel 3', 'Wheel 4' )
	grid

	subplot(2,1,2)
	rpm = d.omega_w;
	plot( c.t, rpm(1,:), 'r', c.t, rpm(2,:), 'g', c.t, rpm(3,:), 'b', c.t, rpm(4,:), 'k' )
	ylabel( 'Wheel rate (rad/sec)' ); xlabel( 'Time (sec)' );
	title( 'Reaction Wheel Angular Velocity' )
	legend('Wheel 1', 'Wheel 2', 'Wheel 3', 'Wheel 4' )
	grid
	zoomrb
	
tile; return
