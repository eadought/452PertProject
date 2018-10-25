function [dudt] = vopODE(t, state, Cd, A, m, JD)
%Variation of Parameters using Gaussian form

    %r_earth = 6378; %km
    mu_earth = 398600; %km^2/s^3
    mu_sun = 132.712e9;
    
    %Initial COES from state
    h0 = state(1);
    ecc0 = state(2);
    inc0 = state(3);
    raan0 = state(4);
    omega0 = state(5);
    theta0 = state(6);
    u0 = omega0 + theta0; %degrees, argument of latitude
    [R0, V0] = COE2VectRad(raan0,inc0,omega0,ecc0,h0,theta0);
    
    %Drag in inertial frame
    H0 = cross(R0, V0);
    r0 = norm(R0);
%     alt = r0-r_earth;
%     W_earth = 72.9211E-6*[0;0;1]; %rad/s, earth angular velocity
%     rho = atmosphere(alt);
%     Vrel = V0 - cross(W_earth,R0); %km/s
%     VrelUnit = Vrel/norm(Vrel); %velocity unit vector
%     P0 = (-.5*rho*1000^2*norm(Vrel)^2*(Cd*A/m)*VrelUnit)/1000; %km/s^2

    %Solar Gravity Acceleration
    JD = JD + t/24/60/60;
    [~,~,Rs] = findSun(JD);
    Rs = Rs';
    Rs_sc = Rs - R0;
    q = dot(R0, (2*Rs - R0))/(norm(Rs)^2);
    Fq = enckeF(q);

    P0 = (mu_sun/norm(Rs_sc)^3)*(Fq*Rs - R0);
    
    %Get LVLH frame in RTN 
    Rvect = R0/norm(R0);
    Nvect = H0/norm(H0);
    temp = cross(Nvect,Rvect);
    Tvect = temp;
    
    %Rotate drag to LVLH
    R = dot(P0,Rvect);
    T = dot(P0,Tvect); 
    N = dot(P0,Nvect);
    
    %%Variate Parameters
    %angular momentum, h
    dhdt = r0*T;
    
    %eccentricity, ecc
    ecc1 = (h0/mu_earth)*sin(theta0)*R;
    ecc2 = (h0^2 + (mu_earth*r0))*cos(theta0);
    ecc3 = mu_earth*ecc0*r0;
    deccdt = ecc1 + ((1/mu_earth/h0)*(ecc2 + ecc3))*T;
    
    %inclination, inc
    dincdt = (r0/h0)*cos(u0)*N;
    
    %right ascension of ascending node, raan
    draandt = (r0*sin(u0)/h0/sin(inc0))*N;
    
    %true anomaly, theta
    theta1 = (h0^2/mu_earth)*cos(theta0)*R;
    theta2 = ((h0^2/mu_earth) + r0)*sin(theta0)*T;
    dthetadt = (h0/r0^2) + ((1/ecc0/h0)*(theta1 - theta2));
    
    %argument of perigee, omega
    dthetapert = (1/ecc0/h0)*(theta1 - theta2);
    domegadt = (-1*dthetapert) - (r0*sin(u0)*N/h0/tan(inc0));
    
    %%state
    dudt = [dhdt;deccdt;dincdt;draandt;domegadt;dthetadt];
end

function [Fq] = enckeF(q)
% Computes the F(q) function used in the encke algorithm
     Fq = ((q^2 - 3*q + 3)/(1 + (1 - q)^(3/2)))*q;
end