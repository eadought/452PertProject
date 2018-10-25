function d2ydt = ode45_doughty ( t ,state )
%Function will perform ODE45 on a 3D vector 
%Input: Vector of initial and final time (T), position (r), velocity (v),
%mu value of orbit (u)
%Output: solution (d2ydt)

u = 398600;
r_mag = norm([state(1) state(2) state(3)]);

%plug in position vector to orbit eq
dxdt = -(u*state(1))/(r_mag^3);
dydt = -(u*state(2))/(r_mag^3);
dzdt = -(u*state(3))/(r_mag^3);

%solution
d2ydt = [ state(4) ; state(5); state(6) ; dxdt ; dydt ; dzdt ];

end