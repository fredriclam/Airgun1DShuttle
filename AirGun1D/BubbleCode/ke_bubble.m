function [dy, dQdt, workrate, dEin] = ke_bubble( ...
    t, y, rho_a, v_a, e_a, p_a, A, physConst, bubbleModel, KErelaxationrate)
    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);
    K = y(5);

    Q       = physConst.Q;
    c_v     = physConst.c_v;
    p_inf   = physConst.p_inf;
    rho_inf = physConst.rho_inf;
    gama    = physConst.gamma;
    c_inf   = physConst.c_inf;

    if strcmpi('quad', bubbleModel)    
        % Surface area and volume factors for hemisphere
        hemisphereFactor = 0.5;
        % Quad bubble rate factor
        rateFactor = 1/4;
    elseif strcmpi('single', bubbleModel)
        hemisphereFactor = 1.0;
        rateFactor = 1;
    end
    
    V = hemisphereFactor * (4/3*pi*R^3);
    Vdot = hemisphereFactor * (4*pi*R^2*Rdot);
    p = E*(gama-1)/V;

    kappa=4000;
    M = 10;
    %M = 25;
    T_inf = 273;
    cv=718;
    Tb = E/(cv*m);
    dQdt = hemisphereFactor*4*pi*R^2*M*kappa*(Tb-T_inf);

    dR = Rdot;
    %dE = A*(e_a + p_a)*v_a - p*Vdot;
    %dE = A*(e_a + p_a)*v_a - p*Vdot - dQdt;
    
    % add turbulent mechanical energy dissipation
    C = 0;
    deltaP = C*rho_inf*abs(Rdot)*Rdot;
    dE = rateFactor*A*(e_a + p_a)*v_a - p*Vdot - dQdt ...
         - hemisphereFactor*4*pi*R^2*Rdot*deltaP;
    workrate = p*Vdot;
    dEin = rateFactor*A*(e_a + p_a)*v_a;
    
    % NEW MECHANICS, replacing above lines
%     KErelaxationrate = 1/(0.0001);
    % Partitioned ratios in choked flow
    kpartition = gama/2;
    epartition = 1/(gama-1);
    partsize = (kpartition + epartition);
    kpartition = kpartition / partsize;
    epartition = epartition / partsize;
    dK = kpartition*dEin - KErelaxationrate*K;
    dE = epartition*dEin + KErelaxationrate*K - p*Vdot - dQdt ...
         - hemisphereFactor*4*pi*R^2*Rdot*deltaP;

    dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;

    %dRdot = 1/R*((p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2);
    %b = 10; 
    %alpha = b*abs(Rdot); %abs(Rdot);
    b = 0;
    alpha=0*0.8; %b*abs(Rdot); %10;
    nu = 0*1e-6;
    % Nonequilibrium factor assuming power-law pressure solution
    a = 0;
    rarefactionFactor = (a + 3) / 3;
    dRdot = 1/R*((rarefactionFactor*p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt ...
        - 3/2*Rdot^2 - alpha*Rdot ...
        - 4 * nu * Rdot / R); % correction from Langhammer and Landro (1996)
    
    dm = rateFactor*A*rho_a*v_a - 0 * 4*pi*R^2 * Rdot * (m/V);

    dy = [dR; dRdot; dm; dE; dK];
end