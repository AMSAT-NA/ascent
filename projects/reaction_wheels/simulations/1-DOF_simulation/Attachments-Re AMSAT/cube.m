% Tetrahedon (see https://en.wikipedia.org/wiki/Tetrahedron)

clear
close all
format short

s2 = sqrt(2);
s3 = sqrt(3);

% Vertices of tetrahedron
%
switch 1
	case 0						% two level edges
		V = [ 1 -1  0  0  0;
			  0  0  1 -1  0;
			[-1 -1  1  1  0] ];
		a = 2;			% edge length

	case 1						% demicube (reaction wheel pyramid)
		V = [ 1  1 -1 -1  0;
			  1 -1  1 -1  0;
			  1 -1 -1  1  0 ];
		a =  2*s2;		% edge length

	case 2						% Vertices on unit sphere 
		V = [2*s2   -s2    -s2    0   0;
			  0    s2*s3 -s2*s3   0   0;
			 -1     -1      -1    3	  0 ]/3;
		a = 2*s2/s3;	% edge length
		
	case 3						% two level edges
		V = [ 1     -1     0     0    0;
			  0      0     1    -1    0;
			[-1/s2 -1/s2  1/s2  1/s2  0] ];
		a = 2;			% edge length


end
 
% mag(V(:,1)-V(:,2)) - a^2				% check edge length

% edge lines
%
E = [ V(:,[1 1 1 2 2 3]); 
	  V(:,[2 3 4 3 4 4]) ];

EC = mean(V(:,1:2),2);					% mid-edge to center
FC = mean(V(:,1:3),2);					% face to center

h = norm(FC) / a;
h = 1/(2*s2*s3);

r = norm(FC-EC) / a;					% radius of cirle on face
r = 1 / (2*s3 );

%h/r - 1/s2								% h/r = 1/sqrt(2)


% coordinates of inscribed circle centered at [0,0,0]'
%	axis thru [0,0,0]' and [1,1,1]'; touches cube 
%   faces at ±[1 -0.5 -0.5], ±[-0.5 1 -0.5], ±[-0.5 -05. 1]
%
ang = pi * [-1:1/3:1]';					% inscribed hexagon
ang = pi * [-1:1/18:1]';				% inscribed circle

circle = s3/s2*[cos(ang), sin(ang), zeros(size(ang)) ]';
circle = [1/s3 0 s2/s3; 0 1 0; -s2/s3 0 1] * circle;
circle1 = [ 1/s2 -1/s2 0;  1/s2  1/s2 0; 0 0 1] * circle;	% 45°
A = [0 -1 0; 1 0 0; 0 0 1];									% 90°
circle2 = A*circle1;
circle3 = A*circle2;
circle4 = A*circle3;


figure
	axis([-1 1 -1 1 -1 1]); hold on
%	plot3( circle1(1,:), circle1(2,:), circle1(3,:), 'k' )
% 	plot3( circle2(1,:), circle2(2,:), circle2(3,:), 'k' )
% 	plot3( circle3(1,:), circle3(2,:), circle3(3,:), 'k' )
 	plot3( circle4(1,:), circle4(2,:), circle4(3,:), 'k' )
	plot3(  V(1,:),  V(2,:),  V(3,:), 'ob' );
	plot3( EC(1), EC(2), EC(3), 'or' );
	plot3( FC(1), FC(2), FC(3), 'og' );
	plot3( E([1 4],:), E([2 5],:), E([3 6],:), 'b' );
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
	axis vis3d
	grid
	rotate3d on


% truncated cone inside cube (axis thru cube center and a vetex)
%
S = 44;					% side of cube containing one wheel

r1 = S/2*s3/s2;			% radius of inscribed circle (centered at cube center)
h1 = r1 * s2;
h = 11.0;
r2 = r1 - h/s2;

A = V(:,1:4)/sqrt(3);		%  body = A * wheel
B = pinv(A);				% wheel = B * body

fprintf( 'S = %4.1f   h = %4.1f   h1 = %4.1f   r1 = %4.1f   r2 = %4.1f\n',	...
		S, h, h1, r1, r2 )
fprintf( 'Angle between flywheel axis and normal of plate = %5.2f°\n', ...
	180/pi*atan( s2 ) );				% face-edge-center angle

	
% 95 x 95 mm mounting plate with 0.25" dia spacers offset 8 mm from edges.
% This is what limits the size of reaction wheel assemby.
%
fprintf( 'center to spacer edge = %5.2f', ( 95/2-8 )*s2 - 0.125*25.4 )
fprintf( '     center to wheel edge = %5.2f\n', r2/s3+(h1+h)*s2/s3 )






