% All variables without subscripts belong to the bubble.
% All variables with subscript _a comes from the airgun.
function dy = bubbleRHS_AllCorrections(y, rho_a, v_a, e_a, p_a, A, physConst, dissp_boolean, beta, Mfac)
    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);

    R_G     = physConst.R_G;
    c_v     = physConst.c_v;
    p_inf   = physConst.p_inf;
    rho_inf = physConst.rho_inf;
    gama    = physConst.gamma;
    c_inf   = physConst.c_inf;

    p = bubblePressure(y, physConst);
    V = 4/3*pi*R^3;
    Vdot = 4*pi*R^2*Rdot;

    kappa=4000;
    M = 10;
    %M = 25;
    T_inf = 273;
    cv=718;
    Tb = E/(cv*m);
    dQdt = 4*pi*R^2*M*kappa*(Tb-T_inf);

    dR = Rdot;
    dE = A*(e_a + p_a)*v_a - p*Vdot - dQdt;
    dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;

    
    if dissp_boolean == true % constant dissipation term
        alpha = beta;
    else
        alpha = beta*abs(Rdot); % dissipation term is proportional to bubble wall velocity
    end
        
    dRdot = (1/R)*(1/rho_inf*(p-p_inf + R*dpdt/c_inf) - 3/2*Rdot^2 - alpha*Rdot);
    
    dm = A*rho_a*v_a*Mfac;

    dy = [dR; dRdot; dm; dE];
end