%% AERO 452 Project 1: Pertubations
% Code by Emily Doughty & Andrew Yazhgur
clear all; close all; clc;

r_earth = 6378;
mu_earth = 398600;
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
%% Object Selection

% Object 1: GEO IntelSat905
% 1 27438U 02027A   18295.82618115 -.00000184  00000-0  00000+0 0  9995
% 2 27438   0.0128 217.0225 0002484 337.8711 109.2763  1.00272477 59949
rp1 =  35776 + r_earth; %km
ra1 =  35797 + r_earth; %km

%Epoch
epoch1 = 18295.82618115;
d1 = 23; y1 = 2018; m1 = 10;
t1 = (epoch1 - floor(epoch1));

% Object 2: LEO Iridium 135
% 1 43070U 17083A   18296.41077455  .00000077  00000-0  20489-4 0  9993
% 2 43070  86.3997  49.5249 0002167  80.4090 279.7350 14.34218808 43715
ra2 = 779 + r_earth; %km
rp2 = 776 + r_earth; %km
%Epoch
epoch2 = 18296.41077455;
d2 = 23; y2 = 2018 ; m2 = 10;
t2 = (epoch1 - floor(epoch2));
%Bstar value from TLE, now this value incorporates more pert affects than
%just drag
bstar2 = 20489e-4; %earthradi^-1

%creat vectors of dates
m = [m1,m2];
d = [d1,d2];
y = [y1,y2];
t = [t1,t2];

%COEs
period = [1436.09,100.4]*60;
ecc = [.0002484,.0002167]; %eccentricity
inc = [0.0128,86.3997]; %inclination (degs)
RAAN = [217.0225,49.5249] ; %right assention of ascending node (deg)
w = [337.8711,80.4090]; %arg of perigee (degs)
MA = [109.2763,279.7350]; %mean anomaly (degs)
n = [1.00272477,14.34218808]; %mean motion (rev/day)

for ii = 1:2
    JD(ii) = UT2JD(y(ii),m(ii),d(ii),t(ii)); 
    a(:,ii) = (mu_earth^(1/3))/(((2*pi*n(ii))/86400)^(2/3)); %semi major axis
	h(:,ii) = sqrt(a(ii)*mu_earth*(1-ecc(ii)^2)); %angular momentum
    %find position and velocity vectors in geocentric frame
    [ R(:,ii) , V(:,ii) , theta(ii) ] = TLE2RV ( MA(ii) , ecc(ii) , h(ii) , RAAN(ii) , inc(ii), w(ii) );
end

%Propagate forward to JD: 18297.0
JD0 = 18297.0;
tspan1 = [0 1+(1-t1)*24*60*60];
tspan2 = [0 (1-t2)*24*60*60];
[~,state_prop1] = ode45( @ode45_doughty , tspan1 , [R(:,1)' V(:,1)'], opt ); %geo
[~,state_prop2] = ode45( @ode45_doughty , tspan2 , [R(:,2)' V(:,2)'], opt ); %leo

Rgeo = [state_prop1(end,1), state_prop1(end,2), state_prop1(end,3)];
Rleo = [state_prop2(end,1), state_prop2(end,2), state_prop2(end,3)];

Vgeo = [state_prop1(end,4), state_prop1(end,5), state_prop1(end,6)];
Vleo = [state_prop2(end,4), state_prop2(end,5), state_prop2(end,6)];

days = 100*24*60*60; %x days in seconds
tspan = [0 days];
stateGeo = [Rgeo Vgeo];
stateLeo = [Rleo Vleo];

[tGeo, dragGeo] = ode45(@CowellODE, tspan, stateGeo, opt, 2.2, 1, 100, JD0, 0);
[tLeo, dragLeo] = ode45(@CowellODE, tspan, stateLeo, opt, 2.2, 1, 100, JD0, 1);

RdragGeo = [dragGeo(:,1), dragGeo(:,2), dragGeo(:,3)];
VdragGeo = [dragGeo(:,4), dragGeo(:,5), dragGeo(:,6)];
RdragLeo = [dragLeo(:,1), dragLeo(:,2), dragLeo(:,3)];
VdragLeo = [dragLeo(:,4), dragLeo(:,5), dragLeo(:,6)];

for i = 1:length(tLeo)
    rdragLeo(i) = norm(RdragLeo(i,:));
end

j = 1;
k = 1;
for i = 2:length(tLeo) - 1
    if(rdragLeo(i) > rdragLeo(i-1) && rdragLeo(i) > rdragLeo(i+1))
        rApoLeo(j) = rdragLeo(i) - r_earth;
        tApoLeo(j) = tLeo(i);
        j = j+1;
    elseif(rdragLeo(i) < rdragLeo(i-1) && rdragLeo(i) < rdragLeo(i+1))
        rPerLeo(k) = rdragLeo(i) - r_earth;
        tPerLeo(k) = tLeo(i);
        k = k+1;
    end
end

plot(tApoLeo/24/60/60, rApoLeo)
hold on
plot(tPerLeo/24/60/60, rPerLeo)
ylabel('Altitude (km)')
xlabel('time (days)')
legend('Apogee', 'Perigee')
title('leo')
    


% [t, Jgeo] = ode45(@OblateCowellODE, tspan, stateGeo, opt, 1,1,1,1,1);
% [t, Jleo] = ode45(@OblateCowellODE, tspan, stateLeo, opt, 1,1,1,1,1);
% 
% plot3(Jleo(:,1), Jleo(:,2), Jleo(:,3))
% hold on
% plot3(Jgeo(:,1), Jgeo(:,2), Jgeo(:,3))



%% Part 1: : Propagate forward two objects (one LEO (<2000 km in semimajor 
% axis and one not in LEO) including Drag, Non-spherical earth, n-body, 
% and solar radiation pressure perturbations when appropriate. 
% Determine the delta-v necessary to correct the perturbations.
% Simplifications
    % Exponential Model for Drag
    % Only J2-J6
    % Only the sun or the moon as extra bodies
    % Same area used for both the SRP and Drag calculations
    
%% Part 2: Pick one of the objects and show that you can correct for the 
% perturbation and continue to correct for a period of time. 
% The tolerances will be variable.

%% Part 3: Pick a new object that uses the perturbations to its advantage 
% to show that some of the delta-v used in part 2 to correct for issues 
% in part 1 are not necessary.
