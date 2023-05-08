clear; clc; close all;

% t = theta for all appearances

% p is our "parameter matrix" containing all given constants and
% information. All of our functions call these constants from the p matrix,
% not necessarily by their names. (EX: the Volume function calls for
% "p(14)" instead of calling for a variable "E")
p = [1.4, .287, 8.4, 101.3, .12, .08, 50*100^(-3), .09, 300, ...
    3*pi/2, pi, 2.8e3, 300, NaN, NaN, 0, 50, 0];
p(14) = (p(6)/(2*p(5))); %calculating E
p(15) = p(4)*volume(pi,p)/p(2)/p(9); %calculating the initial mass

% gam = p(1); %specific heat ratio, unitless
% R= p(2); %gas constant, kJ/kg/K
% r= p(3); %compression ratio, unitless
% P1= p(4); %pressure for t = pi, Kpa
% l= p(5); %connecting rod length, m
% S= p(6); %stroke length, m
% V0= p(7); %volume for t = 0, m^3
% b = p(8); %piston diameter, m
% T1= p(9); %temp for t = pi, K
% thetas= p(10); %t for beginning of heat addition, rad
% thetab= p(11); %angular interval of heat addition, rad
% qin= p(12); %total heat addition, kJ/kg
% Tw= p(13); %cylinder wall temperature, K
% E= p(14); % (s/(2l)), unitless
% M0 = p(15); % mass for t = pi, kg
% C = p(16); % mass blowby constant, unitless
% omega = p(17); %engine speed, rad/s
% hbar = p(18); %convective heat transfer constant, kw/(m^2*K)



% we take symbolic expressions for V and x, and differentiate with respect
% to theta (t) to obtain functions for dV/dt and dx/dt
syms t
symbolicV = (p(7)*(1+(p(3)-1)/(2*p(14))*(1+p(14)*(1-cos(t))- ...
    sqrt(1-p(14)^2*(sin(t))^2))));
dV = diff(symbolicV,t);
dVdt = matlabFunction(dV);

symbolicX = (1/2)*(1-cos(pi*(t-((3*pi)/2))/pi));
dX = diff (symbolicX, t);
dXdt = matlabFunction(dX);


volume(2*pi,p)

%initial guess for ODE45 with initial pressure 101.3 kPa, initial work out
%0 kJ
yInitial = [p(4); 0];

opt = odeset('MaxStep',0.1);
sol = ode45(@(t,y) steamy(t,y,p,dVdt,dXdt), [pi 3*pi], yInitial,opt);

% ---------------------------------------------------------------

for i = 1:length(sol.x)
    vol(i) = volume(sol.x(i), p);
end

% IG equation with new Pressure and volume values. 
% -> To store Temperature Values for plotting. 
% --------------------------------------------------------------
 for i = 1: length(sol.x)
     T(i) = (sol.y(1,i)*volume(sol.x(i),p))/(mass(sol.x(i),p)*(p(2)));
 end 

 % Make 1 by 2 tiled chart layout. 
 % (1,1)-> P-Theta diagram
 % (1,2)-> V-Theta diagram
 % (2,1)-> T-V diagram 
 % (2,2)-> W-theta diagram
 
 tiledlayout(2,2)
% -------------------------------------------------------------- 
% Plot Pressure vs. Crank Angle
 splot = nexttile; 
 plot(sol.x(:), sol.y(1,:))
 hold on 
 title  ("Pressure (Kpa) vs. Crank Angle (Rad)")
 xlabel ("Crank Angle (Rad)")
 ylabel ("Pressure (Kpa)")

% -------------------------------------------------------------- 
% Plot Volume vs. Crank Angle 
 splot = nexttile; 
 plot (sol.x(:),vol(:))
 hold on 
 title  ("Volume Vs. Crank Angle")
 xlabel ("Crank Angle (Rad)")
 ylabel (" Volume (m^3)")
 
% -------------------------------------------------------------- 
% Plot Temperature vs. Crank Angle 
 splot = nexttile;
   plot (sol.x(:),T(:))
   hold on 
   title  ("Temperature (K) vs. Crank Angle (Rad)")
   xlabel ("Crank Angle (Rad)")
   ylabel ("Temperature (K) ")
% -------------------------------------------------------------- 
%   Plot Work vs. Crank Angle 
 splot = nexttile;
   plot(sol.x(:),sol.y(2,:))
   hold on
   title  ("Work (KJ) vs. Crank Angle)")
   xlabel ("Crank Angle (Rad)")
   ylabel ("Work (KJ)")
 
%%
function dy = steamy(t,y,p,dVFunc,dXFunc)
% y(1) = pressure, y(2) = work,out;
% takes input column matrix y = [pressure; workOut; ...] and outputs
% the column vector dy/dt = [dP/dt; dW/dt; d.../dt]

temp = y(1)*volume(t,p)/mass(t,p)/p(2);

dy = zeros(2,1);

dy(1) = -p(1)*(y(1,1)/volume(t,p))*dVFunc(t) + ... 
    (p(1)-1)*mass(t,p)*p(12)*dheatfraction(t,p,dXFunc)/volume(t,p) ...
    -p(1)*p(16)/p(17)*y(1) - (p(1)-1)*p(18)*4/p(8)*(temp-p(13))/p(17);
dy(2) = y(1,1) * dVFunc(t);


end

% --------------------------------------------------------------

function m = mass(t,p)
m = p(15)*exp(-(p(16)/p(17))*(t-pi));
end

% ----------------------------------------------------------------

function V = volume(t,p)
V=(p(7)*(1 + (p(3)-1)/(2*p(14)) * (1+p(14)*(1-cos(t))-sqrt(1-p(14)^2*(sin(t))^2))));
end

% ----------------------------------------------------------------

function x = heatfraction(t,p)

if (pi<=t) && (t<(p(10)))
    x = 0;
elseif ((p(10))<=t) && (t<= ((p(10))+p(11)))
    x = (1/2)*(1-cos(pi*(t-(p(10)))/p(11)));
elseif (( p(10) + p(11) <t)) && (t< 3*pi)
    x = 1;
end

end

% ----------------------------------------------------------------

function dxdt = dheatfraction(t,p,dXFunc)

if (pi<=t) && (t<(p(10)))
    dxdt = 0;
elseif ((p(10))<=t) && (t<= ((p(10))+p(11)))
    dxdt = dXFunc(t);
elseif (( p(10) + p(11) <t)) && (t<= 3*pi)
    dxdt = 0;
end

end