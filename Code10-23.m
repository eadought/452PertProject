%% AERO 452 Project 1: Pertubations
% Code by Emily Doughty & Andrew Yazhgur
clear all; close all; clc;

%% Object Selection

1 27438U 02027A   18295.82618115 -.00000184  00000-0  00000+0 0  9995
2 27438   0.0128 217.0225 0002484 337.8711 109.2763  1.00272477 59949

mu_earth = 398600;
r_earth = 6378;
% Object 1: LEO 
rp1 =  + r_earth; %km
ra1 =  + r_earth; %km
%Epoch
epoch1 = 17331.86035959;
d1 = ; y1 = ; m1 = ;
t1 = (epoch1 - floor(epoch1));

% Object 2: not LEO
ra2 =  + r_earth; %km
rp2 =  + r_earth; %km
%Epoch
epoch2 = 17330.88046138;
d2 = ; y2 = ; m2 = ;
t2 = (epoch1 - floor(epoch2));

%creat vectors of dates
m = [m1,m2];
d = [d1,d2];
y = [y1,y2];
t = [t1,t2];

%COEs
period = [ , ]*60; %seconds
ecc = [ , ]; %eccentricity
inc = [ , ]; %inclination (degs)
RAAN = [ , ] ; %right assention of ascending node (deg)
w = [ , ]; %arg of perigee (degs)
MA = [ , ]; %mean anomaly (degs)
n = [ , ]; %mean motion (rev/day)
for ii = 1:2
    JD(j) = UT2JD(y(j),m(j),d(j),t(j)); 
    a(:,ii) = (mu_earth^(1/3))/(((2*pi*n(ii))/86400)^(2/3)); %semi major axis
	h(:,ii) = sqrt(a(ii)*mu_earth*(1-ecc(ii)^2)); %angular momentum
    %find position and velocity vectors in geocentric frame
    [ R(:,ii) , V(:,ii) , theta(ii) ] = TLE2RV ( MA(ii) , ecc(ii) , h(ii) , RAAN(ii) , inc(ii), w(ii) );
end

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
