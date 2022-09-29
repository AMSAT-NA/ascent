% Tetrahedon (see https://en.wikipedia.org/wiki/Tetrahedron)

clear
close all

rpm = 10000;					% max RPM

rho_br = 8500;					% brass rolled & drawn (kg/m^3)
rho_al = 2720;					% Aluminum 6061 (kg/m^3)
rho_st = 7840;					% Steel 4330 (kg/m^3)
rho_ss = 7800;					% Stainless steel 416 (kg/m^3)

% truncated cone inside cube (axis thru cube center and a vetex)
%		r1		-- radius of cone base
%		r2		-- radius of truncated cone top
%		h1/r1	-- slope of plane cone = sqrt(2)
%		r3,h3	-- hole for motor rotor 
%		r4,h4	-- hole for motor connector clearance (shaft side)
%		r5,h5	-- hole for motor rotor back stop
%
%	           __________________
% 			 /\         |      |
% 			/  \        |      |
% 		   /    \       h2     |
% 		  /      \      |      h1
% 		 /---r2---\  ------    |
% 		/          \    h      |
%	   /-----r1-----\ -----------


% one reaction wheel fits inside an S x S x S cube
%
S = 44e-3;				% cube edge (m)

r1 = S/2*sqrt(3/2);		% radius of inscribed circle (centered at cube center)
h1 = r1 * sqrt(2);

h = 11.0e-3;			% wheel thickness @ widest = 2*h
h2 = h1 - h;
r2 = h2/sqrt(2);

r3 = 10e-3;			h3 =  9e-3;
r4 = 15e-3;			h4 =  h-2e-3;
r5 =  8e-3;			h5 =  h-7e-3;

% h3 + h4 + h5 - 2*h	% should be zero
 

flywheel = {	'cone', r1,  2*h1;
				'cone', r2, -2*h2;
				'cyl',  r3, -h3;
				'cyl',  r4, -h4;
				'cyl',  r5, -h5	};
			
[M,I,V] = rot_inertia( rho_ss, flywheel );	% inertia of flywheel

M_mtr = 25.5e-3;							% Kg = 25.5 g
I_mtr = 3.3e-7;								% Kg m^2 = 3.3 g*cm^2

M = M + M_mtr;								% 1 wheel + 1 motor
I = I + I_mtr;

M_mmt = rho_al * pi*10e-3^2 * 8.5e-3;		% 1 motor mount
M_pyr = rho_al * S^2 * 5e-3;				% 1 pyramid (WAG)
M_plt = rho_al * 1/8*25.4e-3 * ...			% 1/8" mounting plate mass (Kg)
		(95e-3^2 - 2*8e-3^2);		
M_tot = 4*M + 4*M_mmt + M_pyr + M_plt;

% Angular momentum at max RPM
%	(note: N = kg * m / sec^2)
%
L = I * pi/30 *rpm;							% kg-m^2/s or N-m-sec


fprintf( 'S=%5.2f   h=%5.2f   r1=%5.2f   r2=%5.2f   h1=%5.2f', ...
	S*1e3, h*1e3, r1*1e3, r2*1e3, h1*1e3 );

fprintf('   M=%3.0f g   I=%5.0f g-cm^2   Iw=%5.1f mNms   Mt=%3.0f g\n', ...
		M*1e3, I*1e7, L*1e3, M_tot*1e3 );


% Make sure wheels clear the 0.25" dia spacers
% 95 x 95 mm mounting plate with spacers offset 8 mm from edges
%
s2 = sqrt(2);
s3 = sqrt(3);
fprintf( 'center to spacer edge = %5.2f\n', ...
		( 95/2-8 )*s2 - 0.125*25.4 )
fprintf( 'center to wheel edge  = %5.2f\n', ...
		( r2/s3+(h1+h)*s2/s3 )*1e3 ) 
	

% CD disk
%
[M,I,V] = rot_inertia( 1, {'cyl', 60e-3, 1.26e-3} );
M = 16.34e-3;											% (Kg)
rho = M / V;											% (Kg/m^3)
I * rho * 1e7;											% (g-cm^2)




% [M,I,V] = rot_inertia( rho, body )
%	M	 -- mass (Kg)
%	I	 -- rotional inertia (kg-m^2)
%	V	 -- volume (m^3)
%	rho  -- density (Kg/m^3)
%	body -- {'cyl',	  dia (m), h (m) } or
%			{'cone', base (m), h (m) }
%			{'cube', base (m), h (m) }
%
function [M,I,V] = rot_inertia( rho, body )
M=0; I=0; V=0;
for i = 1:size(body,1)
	r = body{i,2};
	h = body{i,3};
	switch body{i,1}
		case 'cyl'
			v = pi * r^2 * h;				% m^3
			m = v * rho;					% kg
			j = 0.5 * m * r^2;				% kg-m^2
		case 'cone'
			v = pi * r^2 * h / 3;			% m^3
			m = v * rho;					% kg
			j = 0.3 * m * r^2;				% kg-m^2
		case 'cube'
			v = r^2 * h;					% m^3
			m = v * rho;					% kg
			j = m * r^2 / 6;				% kg-m^2
		otherwise
			error('Bad body type')
	end
	V = V + v;
	M = M + m;
	I = I + j;
end
return
end

