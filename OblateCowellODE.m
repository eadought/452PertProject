function [dyydt] = OblateCowellODE(t, state, j2, j3, j4, j5, j6)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    R_earth = 6378;
    R = state(1:3);
    r = norm(R);
    mu = 398600;
    J2 = 0.0010826269;
    J3 = -J2*(2.33936E-3);
    J4 = -J2*(1.49601e-3);
    J5 = -J2*(0.20995e-3);
    J6 = J2*0.49941e-3;
    ai = 0;
    aj = 0;
    ak = 0;
    if(j2 == 1)
        ai2 = (-3*J2*mu*(R_earth^2)*R(1)/2/(r^5))*(1-(5*R(3)^2/r^2));
        aj2 = (-3*J2*mu*(R_earth^2)*R(2)/2/(r^5))*(1-(5*R(3)^2/r^2));
        ak2 = (-3*J2*mu*(R_earth^2)*R(3)/2/(r^5))*(3-(5*R(3)^2/r^2));
    
        ai = ai + ai2;
        aj = aj + aj2;
        ak = ak + ak2;
    end
    
    if(j3 == 1)
        ai3 = (-5*J3*mu*(R_earth^3)*R(1)/2/(r^7))*(3*R(3) - (7*R(3)^3/r^2));
        aj3 = (-5*J3*mu*(R_earth^3)*R(2)/2/(r^7))*(3*R(3) - (7*R(3)^3/r^2));
        ak3 = (-5*J3*mu*(R_earth^3)/2/(r^7))*(3*R(3) - (6*R(3)^2 - (7*R(3)^4/r^2) - (3*r^2/5)));
        
        ai = ai + ai3;
        aj = aj + aj3;
        ak = ak + ak3;
    end
    
    if(j4 == 1)
        ai4 = (15*J4*mu*(R_earth^4)*R(1)/(8*r^7))*(1 - (14*R(3)^2/r^2) + (21*R(3)^4/r^4));
        aj4 = (15*J4*mu*(R_earth^4)*R(2)/(8*r^7))*(1 - (14*R(3)^2/r^2) + (21*R(3)^4/r^4));
        ak4 = (15*J4*mu*(R_earth^4)*R(3)/(8*r^7))*(5 - (70*R(3)^2/(3*r^2)) + (21*R(3)^4/r^4));
        
        ai = ai + ai4;
        aj = aj + aj4;
        ak = ak + ak4;
    end
    
    if(j5 == 1)
        ai5 = (3*J5*mu*R_earth^5*R(1)*R(3)/(8*r^9))*(35 - (210*R(3)^2/r^2) + (231*R(3)^4/r^4));
        aj5 = (3*J5*mu*R_earth^5*R(2)*R(3)/(8*r^9))*(35 - (210*R(3)^2/r^2) + (231*R(3)^4/r^4));
        ak5 = ((3*J5*mu*R_earth^5*R(3)/(8*r^9))*(105 - (315*R(3)^2/r^2) + (231*R(3)^4/r^4))) - (15*J5*mu*(R_earth^5)/(8*r^7));
        
        ai = ai + ai5;
        aj = aj + aj5;
        ak = ak + ak5;
    end
    
    if(j6 == 1)
        ai6 = -(J6*mu*R_earth^6*R(1)/(16*r^9))*(35 - (945*R(3)^2/r^2) + (3465*R(3)^4/r^4) - (3003*R(3)^6/r^6));
        aj6 = -(J6*mu*R_earth^6*R(2)/(16*r^9))*(35 - (945*R(3)^2/r^2) + (3465*R(3)^4/r^4) - (3003*R(3)^6/r^6));
        ak6 = -(J6*mu*R_earth^6*R(3)/(16*r^9))*(245 - (2205*R(3)^2/r^2) + (4851*R(3)^4/r^4) - (3003*R(3)^6/r^6));
    
        ai = ai + ai6;
        aj = aj + aj6;
        ak = ak + ak6;
    end
    %propogate orbit, km/s^2
    dxdt = -mu*state(1)/r^3 + ai;
    dydt = -mu*state(2)/r^3 + aj;
    dzdt = -mu*state(3)/r^3 + ak;
    
    dyydt = [state(4); state(5); state(6); dxdt; dydt; dzdt];
    
end

