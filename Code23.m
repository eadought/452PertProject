%% AERO 452 Project 1: Pertubations
% Code by Emily Doughty & Andrew Yazhgur
clear all; close all; clc;

r_earth = 6378;
mu_earth = 398600;
opt = odeset('RelTol',1e-8,'AbsTol',1e-8);
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
%Sat info
%assume mass is average of launch and dry mass of sat (dont know how muh
%fuel is left atm)
m1 = (4723+1984)/2; %kg
A1 = 19.84; %m^2

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
%Sat info
%only info given on mass is launch mass so this is larger than what the
%real mass is atm
m2 = 689; %kg
A2 = 6.90; %m^2

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
% tspan1 = [0 1+(1-t1)*24*60*60];
tspan2 = [0 (1-t2)*24*60*60];
% [~,state_prop1] = ode45( @ode45_doughty , tspan1 , [R(:,1)' V(:,1)'], opt );
[~,state_prop2] = ode45( @ode45_doughty , tspan2 , [R(:,2)' V(:,2)'], opt );

%% Part 1: : Propagate forward two objects (one LEO (<2000 km in semimajor 
% axis and one not in LEO) including Drag, Non-spherical earth, n-body, 
% and solar radiation pressure perturbations when appropriate. 
% Determine the delta-v necessary to correct the perturbations.
% Simplifications
    % Exponential Model for Drag
    % Only J2-J6
    % Only the sun or the moon as extra bodies
    % Same area used for both the SRP and Drag calculations
    
JD = 18297.0;
t = 100; %days 
opt = odeset('RelTol',1e-10,'AbsTol',1e-10,'Event',@terminate);
tspan = [0 24*60*60*t];

% [tnew1_test,state1new_test] = ode45( @ode45_doughty, tspan , [state_prop1(end,1:3)' state_prop1(end,4:6)'] , opt);
[tnew2_test,state2new_test] = ode45( @ode45_doughty, tspan , [state_prop2(end,1:3)' state_prop2(end,4:6)'] , opt);

% [tnew1,state1new] = ode45( @PertProp, tspan , [state_prop1(end,1:3) state_prop1(end,4:6)] , opt , JD);
[tnew2,state2new] = ode45( @PertProp, tspan , [state_prop2(end,1:3) state_prop2(end,4:6)] , opt , JD);

%Check wih plots
% figure
% plot3(state1new_test(:,1),state1new_test(:,2),state1new_test(:,3))
% hold on
% plot3(state1new(:,1),state1new(:,2),state1new(:,3),'*')
% title('GEO')
figure
plot3(state2new_test(:,1),state2new_test(:,2),state2new_test(:,3))
hold on
plot3(state2new(:,1),state2new(:,2),state2new(:,3),'*')
title('LEO')

mue=398600;
% for ii = 1:length(tnew1)
%     [~,e1(ii),inc1(ii),OMEGA1(ii),omega1(ii),~,h1(ii),T(ii)] = RV_COE(state1new(ii,:),398600);
%     rp1(ii)=(h1(ii)^2/mue)*(1/(1+e1(ii)*cosd(0)));
%     ra1(ii)=(h1(ii)^2/mue)*(1/(1+e1(ii)*cosd(180)));
% end
% figure
% subplot(2,2,1)
% hold all
% plot(tnew1/(24*60*60),ra1-6378)
% plot(tnew1/(24*60*60),rp1-6378)
% title('Change in Radius of Apogee & Perigee')
% xlabel('Time (days)'), ylabel('Altitude (km)')
% legend('Apogee','Perigee','Location','Best')
% subplot(2,2,2)
% plot(tnew1/(24*60*60),inc1-inc1(1))
% title('Inclination')
% xlabel('Time (days)'), ylabel('inc (degrees)')
% subplot(2,2,3)
% plot(tnew1/(24*60*60),OMEGA1-OMEGA1(1))
% title('RAAN')
% xlabel('Time (days)'), ylabel('RAAN (degrees)')
% subplot(2,2,4)
% plot(tnew1/(24*60*60),omega1-omega1(1))
% title('Arg of Perigee')
% xlabel('Time (days)'), ylabel('omega (degrees)')
% suptitle('GEO')
for ii = 1:length(tnew2)
    [~,e2(ii),inc2(ii),OMEGA2(ii),omega2(ii),~,h2(ii),T(ii)] = RV_COE(state2new(ii,:),398600);
    rp2(ii)=(h2(ii)^2/mue)*(1/(1+e2(ii)));
    ra2(ii)=(h2(ii)^2/mue)*(1/(1-e2(ii)));
end
figure
subplot(2,2,1)
hold all
plot(tnew2/(24*60*60),ra2-6378)
plot(tnew2/(24*60*60),rp2-6378)
title('Change in Radius of Apogee & Perigee')
xlabel('Time (days)'), ylabel('Altitude (km)')
legend('Apogee','Perigee','Location','Best')
subplot(2,2,2)
plot(tnew2/(24*60*60),inc2-inc2(1))
title('Inclination')
xlabel('Time (days)'), ylabel('inc (degrees)')
subplot(2,2,3)
plot(tnew2/(24*60*60),OMEGA2-OMEGA2(1))
title('RAAN')
xlabel('Time (days)'), ylabel('RAAN (degrees)')
subplot(2,2,4)
plot(tnew2/(24*60*60),omega2-omega2(1))
title('Arg of Perigee')
xlabel('Time (days)'), ylabel('omega (degrees)')
suptitle('LEO')

autoArrangeFigures(2,2,1)

%% Part 2: Pick one of the objects and show that you can correct for the 
% perturbation and continue to correct for a period of time. 
% The tolerances will be variable.

%% Part 3: Pick a new object that uses the perturbations to its advantage 
% to show that some of the delta-v used in part 2 to correct for issues 
% in part 1 are not necessary.
