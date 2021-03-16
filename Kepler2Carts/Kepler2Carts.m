function [X, Y, Z, Vx, Vy, Vz] = Kepler2Carts(alt, ecc, inc, w, nu, RAAN)
%--------------------------------------------------------------------------------------------------------%
%
% 			USAGE: Conversion of Keplerian Classic Orbital Elements into geocentric-equatorial reference system OXYZ
%
% 			AUTHOR: Thameur Chebbi(PhD)		E-MAIL: chebbythamer@gmail.com
%
% 			DATE: 01,Oct,2020
%
% 			DESCIPTION:      This function is created to convert the classical 
%               		     orbital elements to cartesian position and velocity 
%		       		         parameters of any satellite orbit in the geocentric-equatorial
%                            reference system.
%
% 			INPUT:
% 			alt:    Altitude.....................(Km)							
% 			ecc:    Eccentricity											    
% 			inc:	Inclination..................(rad)							
% 			w:	    Argument of perigee..........(rad)	
% 			nu:	    Satellite position...........(rad)							
% 			RAAN:	Right Asc. of Ascending Node.(rad)							
%
% 			OUTPUT:
%			Position Components: 		
% 			[X Y Z]...(Km)
%
%			Velocity Components:
% 			[Vx Vy Vz]...(Km/s) 														
%
%
%%---------------------------- Constants ----------------------------------------------%
mu_earth = 3.986 * 10^5; % Earth Gravitational Constant
re = 6378.1; % earth radius
%
%%--------------------------------------------------------------------------------------
a = alt + re;
p = a*(1-ecc ^2);
r_0 = p / (1 + ecc * cos(nu));
%
%%--------------- Coordinates in the perifocal reference system Oxyz -----------------%
%
% position vector coordinates
x = r_0 * cos(nu);
y = r_0 * sin(nu);
%
%
% velocity vector coordinates
Vx_ = -(mu_earth/p)^(1/2) * sin(nu);
Vy_ = (mu_earth/p)^(1/2) * (ecc + cos(nu));
%
%
%%-------------- the geocentric-equatorial reference system OXYZ ---------------------%
%
% position vector components X, Y, and Z
X = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * x + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * y;
Y = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * x + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * y;
Z = (sin(w) * sin(inc)) * x + (cos(w) * sin(inc)) * y;
% velocity vector components X', Y', and Z'
Vx = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * Vx_ + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * Vy_;
Vy = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * Vx_ + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * Vy_;
Vz = (sin(w) * sin(inc)) * Vx_ + (cos(w) * sin(inc)) * Vy_;
