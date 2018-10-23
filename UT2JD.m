function JD = UT2JD (Y,M,D,t)
%Function will covert given univeral time to julian date
%Inputs: Year (Y), Month (M), Day (D), and time of day fraction (t)
%Outputs: Julian Date (JD)

J0 = 367*Y - floor((7*(Y+floor((M+9)/12)))/4) + floor((275*M)/9)+D+1721013.5;
JD = J0 + (t/24);

end