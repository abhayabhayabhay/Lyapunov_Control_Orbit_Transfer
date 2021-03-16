% Lyapunov Control Orbit Transfer 
% ECH267 - Winter 2021 
% Abhay Negi 

%% README 
%{
    This script calculates...
    How to use: ... 
    Expected simulation time: 
    Notes: 
%}

%% Set Up 

clc
clear
close all 

% Define global variables 
global km m kg s hr Newt rad deg
global aT We Wh
global mu_earth r_earth 
global ef hf

% Define Units 
km       = 1;
m        = km / 1000;
kg       = 1;
s        = 1;
hr       = 3600 * s;
Newt     = kg * m / s^2;
rad      = 1;
deg      = pi * rad / 180;

% Define System Parameters 
mu_earth = 3.986e14 * m^3 * s^(-2);  
r_earth  = 6371 * km;
aT       = 0.010 * m / s^2; % spacecraft acceleration magnitude   

%% Conditions 

% Prompt user to select case to run, see command window when running 
case_select = input('Enter 1 to run altitude change case, Enter 2 to run inclination change case: ') 

% Initial and final conditions
switch case_select 
    case 1
        case_name = 'Altitude Change'
        % Initial Conditions 
        a0       = r_earth + 1000*km;
        e0       = 0;
        inc0     = 0*deg;
        w0       = 0;
        LTAN0    = 0;

        % Final Conditions
        af       = r_earth + 10000*km;
        ef       = 0;
        incf     = 0*deg;
        wf       = 0;
        LTANf    = 0;
    case 2
        case_name = 'Inclination Change'
        % Initial Conditions 
        a0       = r_earth + 1000*km;
        e0       = 0;
        inc0     = 0*deg;
        w0       = 0;
        LTAN0    = 0;

        % Final Conditions
        af       = r_earth + 1000*km;
        ef       = 0;
        incf     = 5*deg;
        wf       = 0;
        LTANf    = 0;
    otherwise 
        error('Invalid case input given.')
end

% Calculate additional initial/final conditions parameters
[h0 e0]  = orbel2vec(mu_earth, a0, e0, inc0, w0, LTAN0);
[hf ef]  = orbel2vec(mu_earth, af, ef, incf, wf, LTANf);

nu       = linspace(0,2*pi,1e3); % time varying parameter for plotting initial orbit
for i = 1:length(nu)
    [x_io(i), y_io(i), z_io(i), u_io(i), v_io(i), w_io(i)] = Kepler2Carts((a0-r_earth)/km, 0, inc0/rad, w0, nu(i), 0);
end
for i = 1:length(nu)
    [x_fo(i), y_fo(i), z_fo(i), u_fo(i), v_fo(i), w_fo(i)] = Kepler2Carts((af-r_earth)/km, 0, -incf/rad, wf, nu(i), 0);
end

r0       = [x_io(1); y_io(1); z_io(1)];
v0       = [u_io(1); v_io(1); w_io(1)];
E0       = norm(v0)^2/2 - mu_earth/norm(r0);
x0       = [r0; v0; h0; e0];

rf       = [x_fo(1); y_fo(1); z_fo(1)];
vf       = [u_fo(1); v_fo(1); w_fo(1)];
Ef       = norm(vf)^2/2 - mu_earth/norm(rf);
xf       = [rf; vf; hf; ef];

%% Simulate 

% note that convergence is sensitive to weights We and Wh 
% if We and Wh are too small, a singularity will occur in the thrust calculation
We       = 1e5*eye(3); 
Wh       = 1e4*eye(3);
t_end    = 120*hr;

tic 
options = odeset('MaxStep',60.0,...
                 'RelTol',1e-4,...
                 'AbsTol',.1);
[t, xs] = ode15s(@xdot, [0, t_end], x0, options); 
[~, T]  = xdot(t.', xs.');
computation_time = toc

%% Post-Process Data 

% Unpack state variables 
x    = xs(:,1);
y    = xs(:,2);
z    = xs(:,3);
u    = xs(:,4);
v    = xs(:,5);
w    = xs(:,6);
hx   = xs(:,7);
hy   = xs(:,8);
hz   = xs(:,9);
ex   = xs(:,10);
ey   = xs(:,11);
ez   = xs(:,12);

% Compute parameters from state variables 
r    = (x.^2 + y.^2 + z.^2).^0.5;
V    = (u.^2 + v.^2 + w.^2).^0.5;
h    = (hx.^2 + hy.^2 + hz.^2).^0.5;
E    = V.^2/2 - mu_earth./r;
e    = (1 + (2*E.*h.^2)/mu_earth^2).^0.5;
evec = [ex, ey, ez];
hvec = [hx, hy, hz];

% Compute Lyapunov function values 
Je   = zeros(length(t),1);
Jh   = zeros(length(t),1);
J    = zeros(length(t),1);
for i = 1:length(t)
    Je(i) = 0.5*(evec(i,:)'-ef)'*We*(evec(i,:)'-ef);
    Jh(i) = 0.5*(hvec(i,:)'-hf)'*Wh*(hvec(i,:)'-hf)/norm(hf)^2;
    J(i)  = Je(i) + Jh(i);      
end

%% State Time Plots 

close all 

save_file_type = '.eps'; % '.png'

% angular momentum plot 
figure(1)
hold on 
plot(t/hr,norm(h0)+0*t,':','linewidth',1)
plot(t/hr,norm(hf)+0*t,'--','linewidth',1)
plot(t/hr,h,'-','linewidth',1)
legend('Initial Ang. Mom.', 'Target Ang. Mom.', 'Trajectory Ang. Mom.','Location', 'Best')
xlabel('Time (hours)')
ylabel('Specific Orbital Angular Momentum, h (km^2/s)')
title([case_name,': Orbital Ang. Mom. Norm. vs. Time'])
box on 
grid on 
saveas(gca,[case_name,' h',save_file_type],'epsc')

% eccentricity plot 
figure(2)
hold on 
plot(t/hr,norm(e0)+0*t,':','linewidth',1)
plot(t/hr,norm(ef)+0*t,'--','linewidth',1)
plot(t/hr,e,'-','linewidth',1)
legend('Initial Eccentricity', 'Target Eccentricity', 'Trajectory Eccentricity','Location', 'Best')
xlabel('Time (hours)')
ylabel('Eccentricity, e (km^2/s)')
title([case_name,': Eccentricity Norm. vs. Time'])
box on 
grid on
saveas(gca,[case_name,' e',save_file_type],'epsc')


% position plot 
figure(3) 
hold on 
plot(t/hr,norm(r0)+0*t,':','linewidth',1)
plot(t/hr,norm(rf)+0*t,'--','linewidth',1)
plot(t/hr,r,'-','linewidth',1)
legend('Initial Position', 'Target Position', 'Trajectory Position','Location', 'Best')
xlabel('Time (hours)')
ylabel('Position (km)')
title([case_name,': Position Norm. vs. Time'])
box on 
grid on
saveas(gca,[case_name,' r',save_file_type],'epsc')


% velocity plot 
figure(4) 
hold on 
plot(t/hr,norm(v0)+0*t,':','linewidth',1)
plot(t/hr,norm(vf)+0*t,'--','linewidth',1)
plot(t/hr,V,'-','linewidth',1)
legend('Initial Velocity', 'Target Velocity', 'Trajectory Velocity','Location', 'Best')
xlabel('Time (hours)')
ylabel('Velocity (km/s)')
title([case_name,': Velocity Norm. vs. Time'])
box on 
grid on
saveas(gca,[case_name, ' v',save_file_type],'epsc')

% energy plot 
figure(5) 
hold on 
plot(t/hr,E0+0*t,':','linewidth',1)
plot(t/hr,Ef+0*t,'--','linewidth',1)
plot(t/hr,E,'-','linewidth',1)
legend('Initial Energy', 'Target Energy', 'Trajectory Energy','Location', 'Best')
xlabel('Time (hours)')
ylabel('Energy (km^2/s^2)')
title([case_name,': Orbital Energy vs. Time'])
box on 
grid on
saveas(gca,[case_name, ' Energy',save_file_type],'epsc')

%% Trajectory Plot

close all

figure('units','normalized','outerposition',[0 0 1 1]) 
hold on
plot3(x,y,z,'linewidth',1)
plot3(x_io,y_io,z_io, '--','linewidth',1)
plot3(x_fo,y_fo,z_fo, '--','linewidth',1)
xlabel('x (km)')
ylabel('y (km)')
xlabel('z (km)')
title([case_name,': Transfer Trajectory'])
legend('Transfer Trajectory', 'Initial Orbit', 'Final Orbit')
axis equal
box on 
grid on

% Earth Trajectory Plot
h = plotearth()
hold on
plot3(x/r_earth,y/r_earth,z/r_earth, 'linewidth', 1, 'color', '#ED1865')
plot3(x_io/r_earth,y_io/r_earth,z_io/r_earth, '--','linewidth',1)
plot3(x_fo/r_earth,y_fo/r_earth,z_fo/r_earth, '--','linewidth',1)
xlabel('x')
ylabel('y')
xlabel('z')
axis equal

%% Supporting Functions 

function [dxdt, T] = xdot(t, x)

    % unpack state vector   
    r = x(1:3,:);
    v = x(4:6,:);
    h = x(7:9,:);
    e = x(10:12,:);
%     mass = x(13,:);
    
    % pull global variables 
    global aT We Wh
    global mu_earth r_earth 
    global Isp g0
    global ef hf
    global xvec 
    xvec = [xvec, x];
    
    Gamma = (1/mu_earth) * skew(v) * skew(r);

    e_term = transpose(e - ef) * We * Gamma;
    h_term = transpose(h - hf) * Wh * skew(r) / norm(hf);
    
    if norm(h-hf) < 100 && norm(e-ef) < 0.001
        u=[0;0;0];
    else
        u = -transpose((e_term + h_term)/norm(e_term + h_term));
    end
    
    %     T = mass*aT*u;
    T = aT*u;    
        
    rdotdot = -mu_earth/norm(r)^3 * r + T;
    hdot = aT * skew(r) * u;
    edot = aT * Gamma * u;
    dxdt = [v; rdotdot; hdot; edot];
%     mdot = -norm(T)/(Isp*g0);
%         dxdt = [v; rdotdot; hdot; edot; mdot];
end

function [h e] = orbel2vec(mu, a, e, i, w, LTAN)
    % mu    = standard gravitational parameter (G*M)  
    % a     = semi-major axis 
    % e     = eccentricity 
    % i     = inclination 
    % w     = argument of periapsis
    % LTAN  = longitude of the ascending node 
    
    h = sqrt(mu*a*(1-e^2)) * transpose(rotmat(w,3)*rotmat(i,1)*rotmat(LTAN,3)) * [0;0;1];
    e = transpose(rotmat(w,3)*rotmat(i,1)*rotmat(LTAN,3)) * [e;0;0];
end

function R = rotmat(angle, axis)
    if axis == 1
        R = [1, 0, 0;
            0, cos(angle), -sin(angle);
            0, sin(angle), cos(angle)];
    elseif axis == 2
        R = [cos(angle), 0, sin(angle);
            0, 1, 0;
            -sin(angle), 0, cos(angle)];
    elseif axis == 3
        R = [cos(angle), -sin(angle), 0;
            sin(angle), cos(angle), 0;
            0, 0, 1];
    else
        error('Error: Axis input must be 1, 2, or 3.')
    end
end

function vskew = skew(v)
    vskew = [    0, -v(3),  v(2)
              v(3),     0, -v(1);
             -v(2),  v(1),    0];
end
% 
% function [X, Y, Z, Vx, Vy, Vz] = Kepler2Carts(alt, ecc, inc, w, nu, RAAN)
% %--------------------------------------------------------------------------------------------------------%
% %
% % 			USAGE: Conversion of Keplerian Classic Orbital Elements into geocentric-equatorial reference system OXYZ
% %
% % 			AUTHOR: Thameur Chebbi(PhD)		E-MAIL: chebbythamer@gmail.com
% %
% % 			DATE: 01,Oct,2020
% %
% % 			DESCIPTION:      This function is created to convert the classical 
% %               		     orbital elements to cartesian position and velocity 
% %		       		         parameters of any satellite orbit in the geocentric-equatorial
% %                            reference system.
% %
% % 			INPUT:
% % 			alt:    Altitude.....................(Km)							
% % 			ecc:    Eccentricity											    
% % 			inc:	Inclination..................(rad)							
% % 			w:	    Argument of perigee..........(rad)	
% % 			nu:	    Satellite position...........(rad)							
% % 			RAAN:	Right Asc. of Ascending Node.(rad)							
% %
% % 			OUTPUT:
% %			Position Components: 		
% % 			[X Y Z]...(Km)
% %
% %			Velocity Components:
% % 			[Vx Vy Vz]...(Km/s) 														
% %
% %
% %%---------------------------- Constants ----------------------------------------------%
% mu_earth = 3.986 * 10^5; % Earth Gravitational Constant
% re = 6378.1; % earth radius
% %
% %%--------------------------------------------------------------------------------------
% a = alt + re;
% p = a*(1-ecc ^2);
% r_0 = p / (1 + ecc * cos(nu));
% %
% %%--------------- Coordinates in the perifocal reference system Oxyz -----------------%
% %
% % position vector coordinates
% x = r_0 * cos(nu);
% y = r_0 * sin(nu);
% %
% %
% % velocity vector coordinates
% Vx_ = -(mu_earth/p)^(1/2) * sin(nu);
% Vy_ = (mu_earth/p)^(1/2) * (ecc + cos(nu));
% %
% %
% %%-------------- the geocentric-equatorial reference system OXYZ ---------------------%
% %
% % position vector components X, Y, and Z
% X = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * x + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * y;
% Y = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * x + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * y;
% Z = (sin(w) * sin(inc)) * x + (cos(w) * sin(inc)) * y;
% % velocity vector components X', Y', and Z'
% Vx = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * Vx_ + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * Vy_;
% Vy = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * Vx_ + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * Vy_;
% Vz = (sin(w) * sin(inc)) * Vx_ + (cos(w) * sin(inc)) * Vy_;
% end
