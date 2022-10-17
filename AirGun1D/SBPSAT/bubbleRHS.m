% All variables without subscripts belong to the bubble.
% All variables with subscript _a comes from the airgun.
function [dy, dQdt, workrate, dEin] = bubbleRHS( ...
    t, y, rho_a, v_a, e_a, p_a, A, physConst, bubbleModel)
    
    Q       = physConst.Q;
    c_v     = physConst.c_v;
    p_inf   = physConst.p_inf;
    rho_inf = physConst.rho_inf;
    gama    = physConst.gamma;
    c_inf   = physConst.c_inf;
    
    %% Leak area
    A_leak = 1e-3; % Approx 1e-3 * pi * 0.30;
    % At t = 0, A is set to approx 50e-3 at t = 0+
%     if t < 0
%         A = A_leak;
%         v_a = sqrt(gama * p_a /rho_a);
%     end
%     if t < 0
%         A = 0;
%     elseif t < 0 || A < A_leak
    if t < 0 || A < A_leak
        A = A_leak;
        v_a = sqrt(gama * p_a /rho_a);
    end


    % Convective heat transfer coefficient
    kappa = 4000;
    % Turbulent magnification factor for heat transfer
    M = 10;
    % Ambient temperature of water
    T_inf = physConst.Tinf;
    % Turbulent mechanical energy dissipation coefficient
    C = 0;
    % Langhammer-Landro dissipation coefficient
    alpha=0.8;
    
    timescaleRelaxation = 0.005; % (unused)
    useWaterFraction = false;
    
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
        if isfield(bubbleModel, 'waterProperties')
            if true
                % Skip
            elseif isfield(bubbleModel.waterProperties, 'fixedWaterFraction')
                useWaterFraction = false;
                c_p = 0.1 * 4e3 + 0.9 * gama*c_v;
                % 10% mass fraction water
                c_v = 0.1 * 4e3 + 0.9 * c_v;
                gama = c_p / c_v;
            else
                useWaterFraction = true;
            end
        end
    else
        error("Unknown bubble model type. Provide 'quad', 'single', " ...
            + "'single-power', or a struct with field 'type' that has " ...
            + "is a string with one of the above values.")
    end
    
    % Check for energy partition model
    useEnergyPartition = false;
    if strcmpi(bubbleModelType, 'partition')
        useEnergyPartition = true;
    end

    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);
    if useEnergyPartition
        K = y(5);
    elseif isfield(bubbleModel, 'waterProperties')
        m_water = y(5);
    end
    
    m_water_rate = 0;
    if useWaterFraction
        % For example, ("psat_T", metadata_reference.discretization.physConst.Tinf-273.15)*1e5
        p_sv = 1.6893e3;
        condensation_coefficient = 0.04;
        m_water_rate = M * sqrt(18.02e-3/(2*pi*(8.314))) ...
            * condensation_coefficient * (p_sv/sqrt(T_inf)) * 4 * pi * R^2;
    else
        m_water_rate = 0;
    end
    
    % Setting for using p data instead of model pressure
    pressureHistoryOverride = false;
    
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
    elseif strcmpi('partition', bubbleModelType)
        % Set all to default:
        hemisphereFactor = 1.0;
        rateFactor = 1;
        rarefactionFactor = 1;
    elseif strcmpi('single-power', bubbleModelType)
        % Single bubble with power-law pressure distribution
        hemisphereFactor = 1.0;
        rateFactor = 1;
        
        % Power p ~ (r/R)^alpha
        pressurePower = -7/3;
        % Surface pressure factor (0 <= this <= 1)
        rarefactionFactor = ...
            (pressurePower * 1 + 3) / 3;
    elseif strcmpi('data-history', bubbleModelType)
        % Set all to default:
        hemisphereFactor = 1.0;
        rateFactor = 1;
        rarefactionFactor = 1;
        pressureHistoryOverride = true;
    else
        error("Unrecognized bubble model type (bubbleModelType)")
    end
    
    if useEnergyPartition
        kpartition = 1;%gama/2;
        epartition = 0;%1/(gama-1);
        partsize = (kpartition + epartition);
        % Compute ratios for sharing incoming enthalpy
        kpartition = kpartition / partsize;
        epartition = epartition / partsize;
    end
    
    c_v_water = 4e3;
    
    % Compute bubble volume
    V = hemisphereFactor * (4/3*pi*R^3);
    Vdot = hemisphereFactor * (4*pi*R^2*Rdot);
    % Compute bubble temperature
    if useWaterFraction
        Tb = E/(c_v*m + c_v_water*m_water);
%         Tb = min([Tb, 20+273]);
    else
        Tb = E/(c_v*m);
    end
    
    if pressureHistoryOverride
        [p, dpdt] = pressureHistory(t);
    elseif useWaterFraction
        % Heat capacity fraction weighting the total energy
%         p = E*(gama-1)/V * c_v*m /(c_v*m + c_v_water*m_water);
        % Ideal gas law pressure assuming incompressible droplets
        massfrac_air = m / (m + m_water);
        p = c_v*(gama-1)*Tb*m/V * massfrac_air;
    else
        % Compute bubble average pressure
        p = E*(gama-1)/V;
    end
    
    deltaP = C*rho_inf*abs(Rdot)*Rdot;
    dQdt = hemisphereFactor*4*pi*R^2*M*kappa*(Tb-T_inf);
    dEin = rateFactor*A*(e_a + p_a)*v_a;
    if useEnergyPartition
        dK = kpartition*dEin - K/timescaleRelaxation;
        dE = epartition*dEin + K/timescaleRelaxation ...
             - p*Vdot - dQdt ...
             - hemisphereFactor*4*pi*R^2*Rdot*deltaP;
    else
        % Water energy influx
        dE_water_in = m_water_rate * c_v_water * T_inf;
        dE = dEin ...
             - p*Vdot - dQdt ...
             - hemisphereFactor*4*pi*R^2*Rdot*deltaP + dE_water_in;
    end
    workrate = p*Vdot;
    
    if pressureHistoryOverride
        % dpdt already retrieved
    elseif useWaterFraction
        % First term in product rule:
        dpdt = (gama-1)*(dE*V-Vdot*E)/V^2 * c_v*m /(c_v*m + c_v_water*m_water);
    else
        dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;
    end

    %dRdot = 1/R*((p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2);
    
    dR = Rdot;
    
    dRdot = 1/R*((rarefactionFactor*p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt ...
        - 3/2*Rdot^2 - alpha*Rdot); % correction from Langhammer and Landro (1996)
    
    dm = rateFactor*A*rho_a*v_a;

    dy = [dR; dRdot; dm; dE];
    
    if useEnergyPartition
        dy = [dy; dK];
    elseif isfield(bubbleModel, 'waterProperties')
        dy = [dy; m_water_rate];
    end
end

% Retrieve pressure from an interpolant on pressure data
function [p, dpdt] = pressureHistory(t)

%% Sum-of-gaussians model
a1 =   5.947e+05;%  (4.891e+05, 7.003e+05)
b1 =    0.003628;%  (0.003542, 0.003714)
c1 =    0.003027;%  (0.002749, 0.003306)
a2 =   5.776e+05;%  (4.878e+05, 6.675e+05)
b2 =    0.007996;%  (0.007006, 0.008987)
c2 =    0.007213;%  (0.006273, 0.008154)
a3 =   8.275e+04;%  (955.5, 1.645e+05)
b3 =      0.2972;%  (0.2874, 0.3069)
c3 =     0.02153;%  (0.01092, 0.03213)
a4 =   2.485e+05;%  (-6238, 5.033e+05)
b4 =      0.0134;%  (-0.000532, 0.02733)
c4 =     0.03032;%  (0.01726, 0.04339)
a5 =   3.838e+17;%  (-3.611e+21, 3.612e+21)
b5 =       8.769;%  (-2778, 2796)
c5 =       1.589;%  (-258.1, 261.3)
a6 =   1.635e+05;%  (-1.797e+05, 5.066e+05)
b6 =     0.02191;%  (-0.1202, 0.164)
c6 =     0.09233;%  (-0.0307, 0.2154)
sum_gaussians = @(x) ...
          a1*exp(-((x-b1)/c1)^2) ...
        + a2*exp(-((x-b2)/c2)^2) ...
        + a3*exp(-((x-b3)/c3)^2) ...
        + a4*exp(-((x-b4)/c4)^2) ...
        + a5*exp(-((x-b5)/c5)^2) ...
        + a6*exp(-((x-b6)/c6)^2);
sum_gaussian_derivatives = @(x) ...
          a1*exp(-((x-b1)/c1)^2)*2/c1^2*(b1-x) ...
        + a2*exp(-((x-b2)/c2)^2)*2/c2^2*(b2-x) ...
        + a3*exp(-((x-b3)/c3)^2)*2/c3^2*(b3-x) ...
        + a4*exp(-((x-b4)/c4)^2)*2/c4^2*(b4-x) ...
        + a5*exp(-((x-b5)/c5)^2)*2/c5^2*(b5-x) ...
        + a6*exp(-((x-b6)/c6)^2)*2/c6^2*(b6-x);

p = arrayfun(sum_gaussians, t);
dpdt = arrayfun(sum_gaussian_derivatives, t);
end