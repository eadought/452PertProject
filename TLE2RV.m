function [R,V, theta] = COE2RV ( MA , ecc, h , RAAN , inc , omega)
% Calculate the position and velocity vector of object given TLE data
% Input: 
%   Mean Anomaly [ MA (degs) }
%   Eccentricity [ ecc ]
%   Angular Momentum [ h (km^2/s) }
%   Right Assention of Ascending Node [ RAAN (degs) }
%   Inclination [ inc (degs) }
%   Argument of Perigee [ omega (degs) }
% Output: 
%   State Vector [ r (km) , v(km/s) ]

u_earth = 398600;
%convert to radians
MA = deg2rad(MA);
RAAN = deg2rad(RAAN);
inc = deg2rad(inc);
omega = deg2rad(omega);
%initial guess for E 
if MA < pi 
    E = MA + ecc;
elseif MA > pi
    E = MA - ecc;
else
end
%iterate to find E
tol = 1e-8; %tolerance
ii = 1; %iteration
ratio = 1;
while abs(ratio) > tol
    fE = E - ecc*sin(E) - MA;
    dfE = 1 - ecc*cos(E);
    ratio = fE/dfE;
    if abs(ratio) > tol
        E = E - ratio;
    else
        E = E;
    end
end
theta = 2*atand(tan(E/2)/(sqrt((1-ecc)/(1+ecc))));
if theta > 360
    theta = theta - 360;
elseif theta < 0
    theta = theta + 360;
end
r_p = cosd(theta); 
r_q = sind(theta);
r_perifocal = (h^2/u_earth)*(1/(1+ecc*cosd(theta)))*[r_p, r_q, 0];
v_p = -sind(theta);
v_q = (ecc + cosd(theta));
v_perifocal = (u_earth/h)*[v_p , v_q , 0];
%convert to geocentric frame
Q_xX = [ ((-sin(RAAN)*cos(inc)*sin(omega))+(cos(RAAN)*cos(omega))), ((-sin(RAAN)*cos(inc)*cos(omega))-(cos(RAAN)*sin(omega))), (sin(RAAN)*sin(inc))  ; ...
         ((cos(RAAN)*cos(inc)*sin(omega))+(sin(RAAN)*cos(omega))),  ((cos(RAAN)*cos(inc)*cos(omega))-(sin(RAAN)*sin(omega))),  (-cos(RAAN)*sin(inc)) ; ...
         (sin(inc)*sin(omega)),                                                 (sin(inc)*cos(omega)),                         (cos(inc))                ];

R = Q_xX*r_perifocal';
V = Q_xX*v_perifocal';
end