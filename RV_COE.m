function [a,e,i,OMEGA,omega,theta,mh,T] = RV_COE(S,mu)
%This function will compute the COE's of an orbit
% mu =  Standard Grav Parameter km^3/s^2
% given a vector s with components of r in the first 3 spots and v in the
% 4-6th spots calculates the COE's

R = [S(1);S(2);S(3)];
V = [S(4);S(5);S(6)];
mV = norm(V); %km/s magnitude of V
mR = norm(R); %km magnitude of R
E = (mV^2)*.5 - mu/mR;

%Semi- Major Axis
a = -1*(mu/(2*E)); %Km

x = mV^2 - mu/mR; 
y = dot(R,V);
ve = (1/mu)*(x*R - y*V); %Eccentricity Vector

% Eccenticitiy
e = norm(ve);
h = cross(R,V); %Angular Momentum
mh = norm(h);
k = [0 0 1]; % K vector

% Inclination
ir = acos(h(3)/ mh); % Inclination (Radians)
i = radtodeg(ir); % Inclination in Degrees


n = cross(k,h); % Ascend Node Vec
mn = norm(n); % Mag of N vector
ih = [1 0 0]; %I vector

% Arguement of Acsending Node
OMEGAR  = acos(dot(ih,n) / mn); % AAN in Radians
OMEGA = radtodeg(OMEGAR); % Arg Of Ascending Node Degree

if n(1,2) <0 
    OMEGA = 360-OMEGA;
end
% Conditional statement for AAN

% Arguement of Periapsis
omegar = acos(dot(n,ve)/(mn*e));%AOP Radians
omega = radtodeg(omegar); %Arg of Periapsis Degree
if ve(3) <0 
    omega = 360 - omega;
end
% conditional statement for AOP

% True Anomaly
thetar = acos(dot(ve,R) / (e*mR)); %True Anomaly Radians
theta = radtodeg(thetar); %True Anomaly in Degrees
if y <0
    theta = 360 - theta;
end
% Conditional statement for True Anomaly
% Calculate Period of the orbit

T = 2*pi*(a^(3/2))/sqrt(mu); %seconds

end
