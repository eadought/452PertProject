function P = dragpert(r_mag,R,v)
%function will calculate the vector of drag pertubation on the s/c

r_alt = r_mag - 6378; %km
rho = atmosphere(r_alt); %kg/m^3
omega_earth = 72.9211e-6*[0 0 1];
cd = 2.2; %coefficent of drag
A = 1; %m^2
m = 100; %kg
V_rel = v - cross(omega_earth,R);
v_rel = norm(V_rel);
uv = V_rel/v_rel;

P = -.5*(cd*A*rho/m)*1000^2*v_rel^2*uv;
P = P/1000;
end


