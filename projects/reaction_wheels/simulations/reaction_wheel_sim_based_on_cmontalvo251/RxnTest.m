%%%Initialize
clear
clc
close all

%%%Globals
global BxI ByI BzI
%%%Setup the IGRF Model
addpath 'igrf/'

%%%Planet Properties
R = 6.371e6; %%meters
M = 5.972e24; %%kg
G = 6.67e-11;  %%%some SI unit
mu = G*M;

%%%Initial Conditions Position and Velocity
altitude = 600*1000; %%meters
x0 = R + altitude;
y0 = 0;
z0 = 0; 
xdot0 = 0;
inclination = 56*pi/180;
semi_major = norm([x0;y0;z0]);
vcircular = sqrt(mu/semi_major);
ydot0 = vcircular*cos(inclination);
zdot0 = vcircular*sin(inclination);
stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];

%%%Intitial Conditions for Attitude and Angular Velocity

phi0 = 0;
theta0 = 0;
psi0 = 0;
%ptp0 = [phi0;theta0;psi0];
%q0123_0 = EulerAngles2Quaternions(ptp0);
p0 = 0.8;
q0 = -0.2;
r0 = 0.3;
%%%Initial conditions of reaction wheels
w10 = 0;
w20 = 0;
w30 = 0;

%state = [x0;y0;z0;xdot0;ydot0;zdot0;q0123_0;p0;q0;r0;w10;w20;w30];


%%%Need time window
period = 2*pi/sqrt(mu)*semi_major^(3/2);
number_of_orbits = 1;
tspan = [0 period*number_of_orbits];
tfinal = period*number_of_orbits;
%tfinal = 200;
timestep = 1.0;
tout = 0:timestep:tfinal;
%stateout = zeros(length(tout),length(state));

%%%This is where we integrate the equations of motion
[tout,stateout] = ode45(@RxnSatTest,tspan,stateinitial);

%%%Loop through stateout to extract Magnetic Field
BxIout = 0*stateout(:,1);
ByIout = BxIout;
BzIout = BxIout;
for idx = 1:length(tout)
    dstatedt = RxnSatTest(tout(idx),stateout(idx,:));
    BxIout(idx) = BxI;
    ByIout(idx) = ByI;
    BzIout(idx) = BzI;
end


%%%Sensor Parameters
lastSensorUpdate = 0;
sensor_params

%%%%Call the Derivatives Routine to initialize vars
%k1 = Satellite(tout(1),state);


%%%Save original State
stateout_original = stateout;



state = [w10;w20;w30];


%%%Globals
global Ir1Bcg  Ir2Bcg  Ir3Bcg n1 n2 n3


%%%RXN WHEELS SIMULATION
%%%Properties
mr = .137;
rr = 43.5/1000;
hr = 24/1000;

%%%One Wheel In each XYZ coordinate
n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];


%%%Inertia Matriz
Idisk = (1/12)*(3*(rr^2)+hr^2);
IrR = mr*[(1/2)*(rr^2) 0 0; 0 Idisk 0; 0 0 Idisk];


%%Transformation from Rxn Wheel frame to body frame of satellite
T1 = Rscrew(n1);
T2 = Rscrew(n2);
T3 = Rscrew(n3);


%%%Inertia of the Rxn Wheel inside the satellite frame
Ir1B = T1'*IrR*T1;
Ir2B = T2'*IrR*T2;
Ir3B = T3'*IrR*T3;


%%%Move the reaction wheels a bit from the cg
r1 = [4;0;0]/1000;
r2 = [0;4;0]/1000;
r3 = [0;0;4]/1000;


%%%Compute the inertia of the reaciotn wheel in the body frame in ref to the cg of the satellite
%%%which means we need to use the parallel axis theorem
sr1 = skew(r1);
Ir1Bcg = Ir1B + mr*(sr1')*sr1;
sr2 = skew(r2);
Ir2Bcg = Ir2B + mr*(sr2')*sr2;
sr3 = skew(r3);
Ir3Bcg = Ir3B + mr*(sr3')*sr3;


%%%Total inertia, I_s is currently 0
Is = 0;
I = Is + Ir1Bcg + Ir2Bcg + Ir3Bcg;

%%%Convert state to kilometers
stateout(:,1:6) = stateout_original(:,1:6)/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);
%q0123out = stateout(:,7:10);
%ptpout = Quaternions2EulerAngles(q0123out);
%pqrout = stateout(:,11:13);
%w123 = stateout(:,14:16);


%{
%%%Total angular momentum
w1 = w123(1);
w2 = w123(2);
w3 = w123(3);
H = Is + Ir1B*w1*n1 + Ir2B*w2*n2 + Ir3B*w3*n3;

%w123dot = [0;0;0];
%}

%%%Make an Earth
[X,Y,Z] = sphere(100);
X = X*R/1000;
Y = Y*R/1000;
Z = Z*R/1000;

%%%Plot 3D orbit
fig = figure();
set(fig,'color','white')
plot3(xout,yout,zout,'b-','LineWidth',4)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on
surf(X,Y,Z,'EdgeColor','none')
axis equal

%%%Plot Magnetic Field
fig2 = figure();
set(fig2,'color','white')
plot(tout,BxIout,'b-','LineWidth',2)
hold on
grid on
plot(tout,ByIout,'g-','LineWidth',2)
plot(tout,BzIout,'r-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Mag Field (nT)')


%%%Plot the angular velocity of the Rxn Wls

fig2 = figure();
set(fig2, 'color', 'white')
plot((tout), w123, 'LineWidth', 2)
grid on
xlabel('Time(sec)')
ylabel('Angular Velocity of Rxn Whls (rad/s)')

toc



