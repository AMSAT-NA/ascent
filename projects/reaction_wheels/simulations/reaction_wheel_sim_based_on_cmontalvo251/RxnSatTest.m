 function dstatedt = Satellite(t,state)
 %%%stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];
 global BxI ByI BzI
x = state(1);
y = state(2);
z = state(3);
%xdot = state(4);
%ydot = state(5);
%zdot = state(6);


%%%Call the magnetic field model
%%%Convert Cartesian x,y,z lattitude and longtitude
phiE = 0;
thetaE = acos(z/rho);
psiE = atan2(y,x);
latitude = 90-thetaE*180/pi;
longitude = psiE*180/pi;
altitude = (rho - R)/1000;
[BxI,ByI,BzI] = igrf('01-Jan-2020', latitude, longtitude, altitude);


%%%inertia parameter
m = 2.6; %%%kilograms

%%%Kinematics
vel = state(4:6);

%%%Gravity Model
planet
r = state(1:3); %% r = [x;y;z]
rho = norm(r);
rhat = r/rho;
Fgrav = -(G*M*m/rho^2)*rhat;

%%%Translational Dynamics
F = Fgrav;
accel = F/m;

dstatedt = [vel;accel];

