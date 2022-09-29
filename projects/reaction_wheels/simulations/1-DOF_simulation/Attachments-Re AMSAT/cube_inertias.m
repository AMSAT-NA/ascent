clear
close all

% Estimate inertias of various Cube sizes
%	MKS is used for all calculations

% see 'Parallel axis theorem' at 
%	https://en.wikipedia.org/wiki/Moment_of_inertia#Parallel_axis_theorem

% 3U estimates from Bob Davis's email:
%
% Ixx ? Iyy ? 36,000 Kg*mm^2 = 36.0e-3 Kg*m^2
% Izz ? 7,200 Kg*mm^2		 =  7.2e-3 Kg*m^2
% M = 4 Kg

%  0U5 dimensions  5 x  5 x  5 cm
%  1U  dimensions 10 x 10 x 10 cm
%  1U5 dimensions 10 x 10 x 15 cm
%  3U  dimensions 10 x 10 x 30 cm
% 12U  dimensions 10 x 20 x 60 cm


% All cubes same density (M*s^3) and CG at their center
%
M = 4/3;								% mass  of 1U (Kg)
s = 0.1;								% sides of 1U (m)


% 1U from homogenious cube
%
I_1U = M * s^2 / 6 * diag( [1 1 1] );	% 1U inertia of a cube (Kg*m^2)
M_1U = M;								% 1U mass (Kg)


% 1.5U from 12 0.5U's					2 x 2 x 3 array of 0U5 cubes
%										0U5's at (±2.5, ±2.5, ±5 ) cm
%												 (±2.5, ±2.5,  0 ) cm
%
M_0U5 = M / 8;
I_0U5 = M_0U5 * (s/2)^2 / 6 * diag( [1 1 1] );		% 0U5 inertia

r = s/4*[ 1  1  2;  1  1 -2;  1 -1 2;  1 -1 -2;
		 -1  1  2; -1  1 -2; -1 -1 2; -1 -1 -2;
		  1  1  0;  1 -1  0; -1  1 0; -1 -1  0 ]';
I_1U5 = 0;
for i=1:size(r,2)
	I_1U5 = I_1U5 +  I_0U5 - M_0U5 * delta(r(:,i))^2;
end
M_0U5 = M_0U5 * size(r,2);


% 3U from 3 1U's						1 x 1 x 3 array of 1U cubes
%										1Us at (0, 0, ±10 ) cm
%										(0, 0,   0 )
r = s*[0 0 1; 0 0 0; 0 0 -1]';	
I_3U = 0;
for i = 1:size(r,2)
	I_3U = I_3U + I_1U - M * delta(r(:,i))^2;
end
M_3U = 3*M;


% 12U from 4 3U's						1 x 2 x 2 array of 3U cubes
%										3Us at (0, ±5, ±15 ) cm
%
J_12U = 0;
r = s/2*[0 1 3; 0 1 -3; 0 -1 3; 0 -1 -3]';
for i = 1:size(r,2)
	J_12U = J_12U + I_3U - M_3U * delta(r(:,i))^2;
end
M_12U = M_3U * size(r,2);


% 12U from 12 1U's						1 x 2 x 6 array of 1U cubes
%										1U's at (0, ±5,  ±5 ) cm
%											    (0, ±5, ±15 ) cm
%											    (0, ±5, ±25 ) cm
%
r = s/2*[0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1;
		 0 1 3; 0 1 -3; 0 -1 3; 0 -1 -3;
		 0 1 5; 0 1 -5; 0 -1 5; 0 -1 -5]';
I_12U = 0;
for i=1:size(r,2)
	I_12U = I_12U + I_1U - M * delta(r(:,i))^2;
end
M_12U = M * size(r,2);


% 12U inertia calculated two ways as sanity check.
% I_12U annd J_12U should be equal



str = '   Ixx=%8.3e   Iyy=%8.3e   Izz=%8.3e\n';
names = {'I_0U5', 'I_1U', 'I_1U5', 'I_3U', 'I_12U'};

for i = 1:length(names)
	fprintf( '\n%3s inertia (Kg*m^2):   ', names{i}(3:end) )
	fprintf( str, diag(eval(names{i})) )
end


% Matrix form of cross product
%	(used to implement 'Parallel axis theorem')
%
function delr = delta( r )
	delr = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
return
end


