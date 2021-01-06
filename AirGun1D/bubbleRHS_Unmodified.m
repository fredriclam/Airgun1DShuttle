% All variables without subscripts belong to the bubble.
% All variables with subscript _a comes from the airgun.
% Reverted Rayleigh-Plesset model with alpha version of Rdot term.
% Same as bubbleRHS (Nov. 17, 2019) but removing the dpdt correction for
% testing analytical asymptotics.
function dy = bubbleRHS_Unmodified(y, rho_a, v_a, e_a, p_a, A, physConst)
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
    %dE = A*(e_a + p_a)*v_a - p*Vdot;
    %dE = A*(e_a + p_a)*v_a - p*Vdot - dQdt;
    
    % add turbulent mechanical energy dissipation
    C = 0;
    deltaP = C*rho_inf*abs(Rdot)*Rdot;
    dE = A*(e_a + p_a)*v_a - p*Vdot - dQdt - 4*pi*R^2*Rdot*deltaP;

    dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;

    %dRdot = 1/R*((p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2);
    %b = 10; 
    %alpha = b*abs(Rdot); %abs(Rdot);
    b = 0;
    alpha=0.8; %b*abs(Rdot); %10;
%     alpha =0;
    dRdot = 1/R*((p-p_inf)/rho_inf + 0*R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2 - alpha*Rdot); % correction from Langhammer and Landro (1996)
    
    dm = A*rho_a*v_a;

    dy = [dR; dRdot; dm; dE];
end