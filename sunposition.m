function Rs = sunposition (JD)
%Function will calculate Rs position vector of sun in km and AU 
%Algorithm 29 from solarlocation notes
%Input: 
%     JD = inital julian date
%Output:
%     Rs_AU = sun position vector in AU
%     Rs_km = sun position vector in km

AU = 149597970.691;
J2000 = 2451545;

Tuti = (JD - J2000)/36525;
lambda_sun = 280.459 + 36000.771*Tuti;
Ttdb = Tuti;
M = 357.5291092 + 35999.0534*Ttdb; 
lambda_elip = lambda_sun + 1.914666471*sind(M) - .019994643*sind(2*M);
rs = 1.000140612 - .016708617*cosd(M) - .000139589*cosd(2*M);
epsilon = 23.439291 - .0130042*Ttdb;

Rs_AU = [cosd(lambda_elip) cosd(epsilon)*sind(lambda_elip) sind(epsilon)*sind(lambda_elip)]*rs;
Rs = Rs_AU*AU;
end