function ic = initialCondtion(rho, v, c, gamma)
    default_arg('c',1);
    default_arg('gamma',1.4);

    p = @(x)rho(x)*c^2/gamma;
    e = @(x)p(x)/(gamma-1)+1/2*rho(x).*v(x).^2;

    ic.rho0 = rho;
    ic.rv0  = @(x)rho(x).*v(x);
    ic.e0   = e;
    ic.p0   = p;
end