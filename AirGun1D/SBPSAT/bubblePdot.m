% All variables without subscripts belong to the bubble.
% All variables with subscript _a comes from the airgun.
function dp = bubblePdot(y, rho_a, v_a, e_a, p_a, A, physConst)
    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);

    gama = physConst.gamma;

    p = bubblePressure(y, physConst);
    Vdot = 4*pi*R^2*Rdot;

    dE = A*(e_a + p_a)*v_a - p*Vdot;

    dp = (gama-1)*dE;
end