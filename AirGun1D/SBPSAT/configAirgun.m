function [physConst, t0, icAirgun, icBubble] = configAirgun(str, ...
    airgunPressure,airgunLength,airgunPortArea,airgunDepth, ...
    airgunCrossSectionalArea, airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength, airgunPortLength, extraOptions)

% Set defaults for last arguments backward compatibility
if nargin == 5
    if strcmpi(str, 'GeneralAirgun')
        error(['Invalid number of input arguments for' ...
               'GeneralAirgun setting.']);
    end
elseif nargin < 11
    error('Invalid number of input arguments.')
elseif nargin < 12
    extraOptions = struct();
end

switch str
    case 'Bolt1500LL'
        % 'realistic' initial conditions that will be used in paper [LW2019].       
        % bubble volume is set equal to the airgun volume 
        % temperature and pressure in the bubble is the same as the ambient properties.
        % The airgun is set to have constant values.
        
        t0 = 0;
        
        p0a_imperial = airgunPressure; % air gun pressure [psi]
        p0a = p0a_imperial * 6894.8; % air gun pressure [Pa]
        physConst.p0a = p0a;
        
        physConst.L = airgunLength; % Length of Airgun in meters

        %cross_sectional_area = 24; % [in^2] cross sectional area is calculated by assuming 1m long (39.3701) for 600in3 airgun.
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
        
        % Cut-off time for fixed opening time for airgun
        physConst.AirgunCutoffTime = 0.100; % time when air gun stops firing [s]
        
                
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
        
        %V_airgun_const = 250 * 1.63871e-5;
        V_airgun_const = 600 * 1.63871e-5;
        icBubble.R = (3/(4*pi) * V_airgun_const)^(1/3);
        %icBubble.R = (3/(4*pi) * V)^(1/3);
        icBubble.Rdot = 0;
        Tvary = T;
        icBubble.m = p*V_airgun_const / (Q*Tvary); 
        icBubble.E = c_v * icBubble.m * Tvary;
        %icBubble.m = p*V_airgun_const / (Q*T); 
        %icBubble.E = c_v * icBubble.m * T;
    
    case 'GeneralAirgun' % Based on the TPS
        % Modified version to accommodate additional parameters relevant to
        % the airgun
        
        % Initial time [s]
        t0 = 0;
        % Airgun pressure conversion
        p0a_imperial = airgunPressure; % Air gun pressure [psi]
        p0a = p0a_imperial * 6894.8;   % Air gun pressure [Pa]
        physConst.p0a = p0a;
        % Airgun length [m]
        physConst.L = airgunLength;

        % Cross-sectional area [in^2]
        crossSectionalAreaImperial = airgunCrossSectionalArea;
        % Compute volume
        airgunLengthImperial = airgunLength/0.0254; % [in]
        V_imperial = crossSectionalAreaImperial * airgunLengthImperial; % [cui]
        physConst.airgunVolume = V_imperial * 1.63871e-5; % Air gun volume [m^3]
        
        %% Water constants
        physConst.rho_inf = 1e3; % density [kg/m^3]
        physConst.pa = 1e5;      % atmospheric pressure [Pa]
        physConst.Tinf = 288;    % temperature (constant) throughout [K]
        depth = airgunDepth;     % depth [m]
        g = 9.8;                 % gravitational acceleration [m/s^2]
        % Compute ambient pressure at depth [Pa]
        physConst.p_inf = physConst.pa + physConst.rho_inf*g*depth;
        
        %% (Air as calorically perfect gas) constants
        physConst.c_v = 718;    % heat capacity of air at constant volume [J/kgK]
        physConst.c_inf = 1482; % speed of sound in water [m/s]
        physConst.Q = 287.06;   % specific gas constant for dry air [J/kgK]
        physConst.gamma = 1.4;  % ratio of heat capacities for dry air
        
        %% Compute airgun initial conditions
        
        if isfield(extraOptions, 'TInitial')
            T = extraOptions.TInitial;
            disp('Using extra option -- TInitial -- for bubble & chamber');
        else
            T = physConst.Tinf;
        end
        if isfield(extraOptions, 'gasConstant')
            physConst.Q = extraOptions.gasConstant;
            disp('Using extra option -- gasConstant Q');
        elseif isfield(extraOptions, 'Q')
            physConst.Q = extraOptions.Q;
            disp('Using extra option -- gasConstant Q');
        end
        if isfield(extraOptions, 'gamma')
            physConst.gamma = extraOptions.gamma;
            disp('Using extra option -- gamma');
        end
        if isfield(extraOptions, 'c_v')
            physConst.c_v = extraOptions.c_v;
            disp('Using extra option -- c_v');
        end
        
        p = physConst.p0a;
        Q = physConst.Q;
        c_v = physConst.c_v;
        rho = p/(Q*T);
        e = rho*c_v*T;
        % Set uniform, quiescent airgun initial conditions
        icAirgun.rho0 = @(x)0*x + rho;
        icAirgun.rv0  = @(x)0*x;
        icAirgun.e0   = @(x)0*x + e;
        icAirgun.p0   = @(x)0*x + p;
        
        if isfield(extraOptions, 'IC')
            try
                IC = extraOptions.IC;
                icAirgun.rho0 = IC.rho0;
                icAirgun.rv0  = IC.rv0;
                icAirgun.e0   = IC.e0;
                icAirgun.p0   = IC.p0;
            catch
                error("Could not unpack initial condition (IC) " ...
                    + "functions rho0, rv0, e0, and p0.")
            end
        end
        
        %% Set bubble initial conditions and gas air parameters
        p = physConst.p_inf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        % Convert to metric the initial bubble volume argument
        % Note: LeightonWatson/airgun1d/HEAD uses a fixed presribed volume
        % 600, off by just a little from the chamber volume [cui]
        bubbleInitialVolumeMetric = bubbleInitialVolume * 1.63871e-5; 
        icBubble.R = (3/(4*pi) * bubbleInitialVolumeMetric)^(1/3);
        icBubble.Rdot = 0;
        icBubble.m = p*bubbleInitialVolumeMetric / (Q*T); 
        icBubble.E = icBubble.m * c_v * T;
        
        is_icBubbleModified = false;
        if isfield(extraOptions, 'R')
            icBubble.R = extraOptions.R;
            disp('Using extra option -- R');
            is_icBubbleModified = true;
        end
        if isfield(extraOptions, 'Rdot')
            icBubble.Rdot = extraOptions.Rdot;
            disp('Using extra option -- Rdot');
            is_icBubbleModified = true;
        end
        if isfield(extraOptions, 'pb')
            pb = extraOptions.pb;
            icBubble.m = pb*(4/3*pi*icBubble.R^3) / (Q*T); 
            icBubble.E = icBubble.m * c_v * T;
            disp('Using extra option -- pb');
            is_icBubbleModified = true;
        end
        
        % Convert cross sectional area from [in^2] to [m^2]
        physConst.crossSectionalArea = airgunCrossSectionalArea*0.00064516;
        % Set left pressure = right pressure
        physConst.p_R0 = p0a; % [Pa]
        
        % Set fixed shuttle assembly mass [kg]
        physConst.shuttleAssemblyMass = 63 * .454;
        if isfield(extraOptions, 'shuttleAssemblyMass')
            physConst.shuttleAssemblyMass = ...
                extraOptions.shuttleAssemblyMass;
            disp('Using extra option -- shuttleAssemblyMass');
        end
        
        % Convert port area from [in^2] to [m^2]
        physConst.APortTotal = airgunPortArea*0.00064516;
        % From measurement:
        physConst.APortTotal = 90 * (0.0254)^2; % 0.0581 m^2
        % Set fixed operating chamber length [m]
        physConst.operatingChamberLength = 87.6e-3;
        
        % Set firing chamber and operating chamber profiles A(x)
        physConst.airgunFiringChamberProfile = ...
            airgunFiringChamberProfile;
        physConst.airgunOperatingChamberProfile = ...
            airgunOperatingChamberProfile;
        % Set cupped piston area
        physConst.shuttle_area_left = 0.04735; % [m^2]
        % Set front-of-shuttle (flange) projected area
        physConst.shuttle_area_right = 0.06155; % [m^2]
        % Rear of shuttle: projected area of piston, minus the shaft area
        % Estimate of shaft area:
        shaft_area = pi/4 * (2.1 * 0.0254)^2;
        physConst.shuttle_area_right_rear = ...
            physConst.shuttle_area_right - shaft_area; % [m^2]
        % Lead-in length (acceleration distance) where shuttle can move
        % without exposing the air
        physConst.portLead = 0.0; % [m]
        
        % (Unused) thickness of the flange in the operating chamber
        physConst.flangeDepth = 3 * 0.0254; % [m]
        % Approximate the flange ID to be equal to the chamber
        
        physConst.plugVolumeFn = @(xi) error("Deprecated function"); ...
            % physConst.crossSectionalArea * (xi + physConst.flangeDepth);

        % Set shuttle penalty parameter
        %   (~1e11 makes sense based on linear elasticity, although
        %    this should typically not activate because of the gas region
        %    cushioning the closing shuttle III.)
        physConst.shuttleBdryPenaltyStrength = shuttleBdryPenaltyStrength;
        
        physConst.airgunPortLength = 84.6e-3/0.0254; %[in]
        physConst.midChamberLength = 3.112 * 0.0254;
        physConst.midChamberArea = physConst.shuttle_area_left;
        physConst.OpRearOrificeArea =  1.37e-6; % (0.052"D)
        
        if isfield(extraOptions, 'dampingConstant')
            physConst.dampingConstant = ...
                extraOptions.dampingConstant;
            disp('Using extra option -- dampingConstant');
        end
        if isfield(extraOptions, 'dampingQuadraticConstant')
            physConst.dampingQuadraticConstant = ...
                extraOptions.dampingQuadraticConstant;
            disp('Using extra option -- dampingQuadraticConstant');
        end
        if isfield(extraOptions, 'OpRearOrificeArea')
            physConst.OpRearOrificeArea = ...
                extraOptions.OpRearOrificeArea;
            disp('Using extra option -- OpRearOrificeArea');
        end
        
    otherwise
        error("in configAirgun: unknown parameters key");
end
end