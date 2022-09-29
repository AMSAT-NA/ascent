clear
close all

rng('default')

% L converts from wheel momentums to a resulting momentum vector, i.e
%	x = L*w (also w = pinv(L)*x)
%  
L = [ 1 -1 -1  1;					% axes of the 4 wheels
	  1  1 -1 -1;
	  1 -1  1 -1 ] / sqrt(3);
  
L(:,4) = 0; 						% disable wheel 4

% generate random points on a sphere
%
N = 50e3;
TH = 2*pi*rand(1,N);
PH = asin(-1+2*rand(1,N));
[X,Y,Z] = sph2cart(TH,PH,1);
x = [X;Y;Z];


% x are momentum vectors we might want to use to rotate the satellite
% or maintain the satellite orientation
%
w = pinv(L)*x;						% convert to wheel momentums

w_max = max(max(abs(w)));			% max wheel momentum
disp( w_max )

r = 400/6000;						% min RPM / max RPM


% don't plot points if wheel RPM < r*w_max
%
i = find( abs(w(1,:)) < r*w_max );	w(:,i) = NaN;
i = find( abs(w(2,:)) < r*w_max );	w(:,i) = NaN;
i = find( abs(w(3,:)) < r*w_max );	w(:,i) = NaN;
% i = find( abs(w(4,:)) < r*w_max );	w(:,i) = NaN;


% percent points in dead bands
%
display( 100*sum(isnan(sum(w)))/N )


% back to momentum vectors 
%	(y==x except where abs(wheel momentum) < r*w_max)
%
y = L*w;			


% check that x==y except where abs(wheel momentum) < r*w_max)
%
[~,j] = find( ~isnan(sum(w)) );
all(all( abs( x(:,j)-y(:,j)) < 10*eps ))


% generate elipses in center of dead bands
%
ang = pi * [-1:1/18:1]';
c = cos(ang);
s = sin(ang);
z = zeros(size(ang));
bandx = 1.3*L*[z c s z]';			% wheel 1 & 4 = 0
bandy = 1.3*L*[c z s z]';			% wheel 2 & 4 = 0
bandz = 1.3*L*[c s z z]';			% wheel 3 & 4 = 0

	
% Axis of a dead band, for example bandx, can be found with:
%
% v1 = L*[0 c s 0]' with c = 1 & s = 0
% v2 = L*[0 c s 0]' with c = 0 & s = 1
% axis of bandx = cross( v1, v2 ) / 2;
%
% axis of bandx = [0  1  1]';
% axis of bandy = [1 -1  0]';
% axis of bandz = [1  0 -1]';

figure
	axis(1.1*[-1 1 -1 1 -1 1]); hold on
	plot3( y(1,:), y(2,:), y(3,:), '.c', 'markersize', .125 )
	plot3( [0  0], [1 -1], [1 -1], '--r', 'LineWidth', 2 )
	plot3( [1 -1], [-1 1], [0  0], '--g', 'LineWidth', 2 )
	plot3( [1 -1], [0  0], [-1 1], '--b', 'LineWidth', 2 )
	plot3( bandx(1,:), bandx(2,:), bandx(3,:), 'r', 'LineWidth', 1  )
  	plot3( bandy(1,:), bandy(2,:), bandy(3,:), 'g', 'LineWidth', 1  )
  	plot3( bandz(1,:), bandz(2,:), bandz(3,:), 'b', 'LineWidth', 1  )
	axis vis3d
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
	grid
	rotate3d on

