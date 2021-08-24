% All variables without subscripts belong to the bubble.
% All variables with subscript _a comes from the airgun.
function [dy, dQdt, workrate, dEin] = bubbleRHS( ...
    t, y, rho_a, v_a, e_a, p_a, A, physConst, bubbleModel)
    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);

    Q       = physConst.Q;
    c_v     = physConst.c_v;
    p_inf   = physConst.p_inf;
    rho_inf = physConst.rho_inf;
    gama    = physConst.gamma;
    c_inf   = physConst.c_inf;

    % Convective heat transfer coefficient
    kappa = 4000;
    % Turbulent magnification factor for heat transfer
    M = 10;
    % Ambient temperature of water
    T_inf = 273 + 15; % Changed
    % Turbulent mechanical energy dissipation coefficient
    C = 0;
    % Langhammer-Landro dissipation coefficient
    alpha=0.8;
    % Power p ~ (r/R)^alpha
    pressurePower = -7/3;
    timescaleRelaxation = 0.1;
    
    % Allow bubble model to be specified as a character-array or as a
    % struct with field "type".
    if isstr(bubbleModel)
        bubbleModelType = bubbleModel;
    elseif isstruct(bubbleModel)
        bubbleModelType = bubbleModel.type;
        if isfield(bubbleModel, 'M')
            M = bubbleModel.M;
        end
        if isfield(bubbleModel, 'alpha')
            alpha = bubbleModel.alpha;
        end
        if isfield(bubbleModel, 'pressurePower')
            pressurePower = bubbleModel.pressurePower;
        end
        if isfield(bubbleModel, 'timescaleRelaxation')
            timescaleRelaxation = bubbleModel.timescaleRelaxation;
        end
    else
        error("Unknown bubble model type. Provide 'quad', 'single', " ...
            + "'single-power', or a struct with field 'type' that has " ...
            + "is a string with one of the above values.")
    end
    
    if strcmpi('quad', bubbleModelType)    
        % Surface area and volume factors for hemisphere
        hemisphereFactor = 0.5;
        % Quad bubble rate factor
        rateFactor = 1/4;
        % Disable rarefaction factor
        rarefactionFactor = 1;
    elseif strcmpi('single', bubbleModelType)
        % Set all to default:
        hemisphereFactor = 1.0;
        rateFactor = 1;
        rarefactionFactor = 1;
    elseif strcmpi('single-power', bubbleModelType)
        % Single bubble with power-law pressure distribution
        hemisphereFactor = 1.0;
        rateFactor = 1;
        
        % Relaxation timescale
%         tEqu = 1/100;
        % Exponential factor, activating after timescale of expansion
%         exponentialFactor = exp(-max(0, t - timescaleRelaxation) / tEqu);
%         tanhFactor = 1/2 + 1/2*tanh(-(t-timescaleRelaxation)/tEqu);
%         gaussianFactor = exp(-((t/timescaleRelaxation)^2));
        % Surface pressure factor (0 <= this <= 1)
        rarefactionFactor = ...
            (pressurePower * 1 + 3) / 3;
    end
    
    % Compute bubble temperature
    Tb = E/(c_v*m);
    % Compute bubble volume
    V = hemisphereFactor * (4/3*pi*R^3);
    Vdot = hemisphereFactor * (4*pi*R^2*Rdot);
    % Compute bubble average pressure
    p = E*(gama-1)/V;

    deltaP = C*rho_inf*abs(Rdot)*Rdot;
    dQdt = hemisphereFactor*4*pi*R^2*M*kappa*(Tb-T_inf);
    dE = rateFactor*A*(e_a + p_a)*v_a - p*Vdot - dQdt ...
         - hemisphereFactor*4*pi*R^2*Rdot*deltaP;
    workrate = p*Vdot;
    dEin = rateFactor*A*(e_a + p_a)*v_a;

    dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;

    %dRdot = 1/R*((p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2);
    
    dR = Rdot;
    
    dRdot = 1/R*((rarefactionFactor*p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt ...
        - 3/2*Rdot^2 - alpha*Rdot); % correction from Langhammer and Landro (1996)
    
    dm = rateFactor*A*rho_a*v_a;

    dy = [dR; dRdot; dm; dE];
end