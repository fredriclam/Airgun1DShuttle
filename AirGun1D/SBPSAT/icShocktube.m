function ic = icShocktube()
    gamma = 1.4;

    rho_l = 1;
    rho_r = 0.125;
    p_l = 1;
    p_r = 0.1;
    c = 1; %speed of sound

    v   = @(x)0*x;
    rho = @(x)(x<=0)*rho_l + (x>0)*rho_r;
    p   = @(x)(x<=0)*p_l + (x>0)*p_r;
    e   = @(x)p(x)/(gamma-1);

    ic.rho0 = rho;
    ic.rv0  = @(x)rho(x).*v(x);
    ic.e0   = e;
    ic.p0   = p;
end