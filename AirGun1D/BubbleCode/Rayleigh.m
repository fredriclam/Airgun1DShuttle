function dYdt = Rayleigh(t, Y, params)
% dYdt = Rayleigh(t, Y, params);

% This function file solves the Herring equation which, similarly to the
% Keller equation, assumes a constant speed of sound in the fluid. This
% file solves the modified Herring equation given by equation (6) in 
% Vokurka, (1986), Comparison of Rayleigh's, Herring's, and Gilmore's
% models of gas bubbles
 
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

%% COMPUTE REQUIRED QUANTITIES %%

V = 4/3*pi*R^3; % bubble volume
dVdt = 4*pi*R^2*U; % rate of change of bubble volume

%% EVOLVE SOLUTION %%

dTdt = 1/(cv*m0) * (-4*pi*R^2*M*kappa*(Tb-T_inf) - pb*dVdt); % bubble temperature

dpdt = (1/V^2) * (m0*Q*dTdt*V - m0*Q*Tb*dVdt); % bubble pressure

dRdt = U; % bubble radius

%dUdt = (1/R) * ((1/rho_infty) * (pb-p_infty+R/c_infty*dpdt) - 3/2*U^2); % bubble wall velocity
dUdt = (1/R) * ((1/rho_infty) * (pb-p_infty) - 3/2*U^2); % bubble wall velocity

dYdt = [dRdt; dUdt; dTdt; dpdt]; % evolve solution