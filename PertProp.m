function d2ydt = PertProp ( t ,state , JD)
%Function will perform ODE45 on a 3D vector including Drag, Non-spherical 
%earth, n-body, and solar radiation pressure
%Input: Vector of initial and final time (T), position (r), velocity (v),
%mu value of orbit (u)
%Output: solution (d2ydt)

u = 398600;
r_mag = norm([state(1) state(2) state(3)]);
R = state(1:3);
V = state(4:6);
%Update JD each step
JD = JD + t/86400;
% SRP
Rs = sunposition(JD);
[srp,~] = srp_funct( R , Rs , 1 );
% DRAG
P = dragpert(r_mag,R,V);
% Total Pertubational Affect
PERT = srp + P;

%plug in position vector to orbit eq
dxdt = -(u*state(1))/(r_mag^3) + PERT(1);
dydt = -(u*state(2))/(r_mag^3) + PERT(2);
dzdt = -(u*state(3))/(r_mag^3) + PERT(3);

%solution
d2ydt = [ state(4) ; state(5); state(6) ; dxdt ; dydt ; dzdt ];

end