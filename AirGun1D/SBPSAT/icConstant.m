function ic = icConstant(p0)
    c = 1;
    gamma = 1.4;

    p = @(x) x*0 + p0;
    rho = @(x)p(x)*gamma/c^2;
    v = @(x)x*0;
    e = @(x)p(x)/(gamma-1)+1/2*rho(x).*v(x).^2;

    ic.rho0 = rho;
    ic.rv0  = @(x)rho(x).*v(x);
    ic.e0   = e;
    ic.p0   = p;
end