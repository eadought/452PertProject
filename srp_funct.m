function [srp,f] = srp_funct( R , Rs , on )
%Function will perform ODE45 on a 3D vector 
%Input: 
%     R = position vector
%     cr = reflectivity
%     Asun = exposed area to sun
%     Psr = solar radiation pressure
%     Rs = distance earth to the sun
%Output:
%     Psr = Solar Radiation Pressure
%shadow funciton

r_earth = 6378;
Psr = 4.57e-6; %N/m^2
cr = 1.2;
Asun = 1; 
r = norm(R);
rsmag = norm(Rs);
rsmag = norm(Rs);
theta = acos(dot(Rs,R)/(rsmag*r));
theta1 = acos(r_earth/r);
theta2 = acos(r_earth/rsmag);
Rs_sun = Rs - R;
rs_sun = norm(Rs_sun);
m = 100;

if on == 1
    if (theta1+theta2 <= theta)
        f = 0;
    else 
        f = 1;
    end
elseif on == 0
    f = 1;
else
end

%calculate solar rad pressure force for each step
srp = -(Psr*cr*Asun/m)*(Rs_sun/rs_sun)*f;
end