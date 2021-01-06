function dYdt = Gilmore(t, Y, params);
% dYdt = gilmore(t, Y, params);

% This function file solves the Gilmore equation. This is different to the
% Keller or Herring equation because the speed of sound varies with the
% pressure, rather than assumed to be fixed. The Gilmore equation is given
% by equations (3), (4), and (5) in Vokurka, (1986), Comparison of 
% Rayleigh's, Herring's, and Gilmore's models of gas bubbles

%% LOAD VARIABLES %%

R = Y(1); %bubble radius [m]
U = Y(2); %bubble velocity [m/s]
Tb = Y(3); % bubble temperature [K]
pb = Y(4); % bubble pressure [Pa]

%% LOAD PARAMETERS %%

rho_infty = params(1); %ambient density
p_infty = params(2); %ambient pressure  
p0 = params(3); %initial pressure of bubble
R0 = params(4); %initial radius of bubble
k = params(5); %polytropic index
c_infty = params(6); %speed of sound in fluid

cv = 718; % heat capacity
M = 10; % coefficient of turbulent surface area
kappa = 4000; % heat capacity coefficient
T_inf = 288; % ambient temperature
Q = 287.06; % gas constant
m0 = params(7); % mass in bubble

B = params(8); %constant from Tait equation of state [Pa]
n = params(9); %constant from Tai equation of state [ ] 

%% COMPUTE REQUIRED QUANTITIES %%

V = 4/3*pi*R^3; % bubble volume
dVdt = 4*pi*R^2*U; % rate of change of bubble volume

dTdt = 1/(cv*m0) * (-4*pi*R^2*M*kappa*(Tb-T_inf) - pb*dVdt); % bubble temperature

dpdt = (1/V^2) * (m0*Q*dTdt*V - m0*Q*Tb*dVdt); % bubble pressure

c = c_infty * (((pb+B)/(p_infty+B))^((n-1)/(2*n))); %speed of sound at bubble wall

H = n/(n-1)*(p_infty+B)/rho_infty * (((pb+B)/(p_infty+B))^((n-1)/n) - 1); %enthalpy
dHdt = dpdt/rho_infty * ((pb+B)/(p_infty+B))^(-1/n); %rate of change of enthalpy

%% EVOLVE SOLUTION %%

dRdt = U;
dUdt = (H*(1+U/c) + R/c*dHdt*(1-U/c) - 3/2*U^2*(1-U/(3*c))) / (R*(1-U/c));

dYdt = [dRdt; dUdt; dTdt; dpdt];
