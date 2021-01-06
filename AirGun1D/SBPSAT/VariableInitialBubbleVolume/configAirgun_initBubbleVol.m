function [physConst, t0, icAirgun, icBubble] = configAirgun_initBubbleVol(str,initBubbleVol,airgunPressure,airgunLength,airgunPortArea,airgunDepth)

switch str
            
    case 'Bolt1500LL'
        % 'realistic' initial conditions that will be used in paper.        
        % bubble volume is set equal to the airgun volume 
        % temperature and pressure in the bubble is the same as the ambient properties.
        % The airgun is set to have constant values.
        
        t0 = 0;
        
        p0a_imperial = airgunPressure; % air gun pressure [psi]
        p0a = p0a_imperial * 6894.8; % air gun pressure [Pa]
        physConst.p0a = p0a;
        
        physConst.L = airgunLength; % Length of Airgun in meters

        cross_sectional_area = airgunPortArea;
        airgunLengthImperial = physConst.L/0.0254; % length in inches
        V_imperial = airgunLengthImperial * cross_sectional_area;
        V = V_imperial * 1.63871e-5; % air gun volume [m^3]
        
        A_imperial = airgunPortArea; % air gun port area [in^2]
        A = A_imperial * 6.4516e-4; % air gun port area [m^2]
        physConst.A = A;
        
        physConst.rho_inf = 1e3; % density [kg/m^3]
        physConst.pa = 1e5; % atmospheric pressure [Pa]
        depth = airgunDepth; % depth [m]
        g = 9.8; % gravitational acceleration [m/s^2]
        physConst.p_inf = physConst.pa + physConst.rho_inf*g*depth; % ambient pressure at depth [Pa]
        
        physConst.c_v = 718; % heat capacity of air at constant volume [J/kgK]
        physConst.c_inf = 1482; % speed of sound in water [m/s]
        physConst.Q = 287.06; % specific gas constant for dry air [J/kgK]
        physConst.R_G = physConst.Q; % DUPLICATE
        physConst.Tinf = 288; % temperature assumed constant throughout the system [K]
        physConst.gamma = 1.4; % ratio of heat capacities for dry air
        
        physConst.AirgunCutoffTime = 0.01; % time when air gun stops firing
        
        % Air gun
        p = physConst.p0a;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        rho = p/(Q*T);
        e = c_v*rho*T;
        
        icAirgun.rho0 = @(x)0*x + rho;
        icAirgun.rv0  = @(x)0*x;
        icAirgun.e0   = @(x)0*x + e;
        icAirgun.p0   = @(x)0*x + p;
        
        % Bubble
        p = physConst.p_inf;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        V_airgun_const = initBubbleVol * 1.63871e-5;
        icBubble.R = (3/(4*pi) * V_airgun_const)^(1/3);
        %icBubble.R = (3/(4*pi) * V)^(1/3);
        icBubble.Rdot = 0;
        icBubble.m = p*V_airgun_const / (Q*T); 
        icBubble.E = c_v * icBubble.m * T;
        
              
    otherwise
        error();
    end
end