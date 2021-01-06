function ic = icDensityBell(v0, rho0, A, x0, d)
    default_arg('v0', 0);
    default_arg('rho0', 0.1);
    default_arg('A', 0.2);
    default_arg('x0', 0);
    default_arg('d', sqrt(0.05));

    rho = @(x)x*0 + rho0 + A*exp(-(x-x0).^2/d^2);
    v = @(x)x*0 + v0;
    ic = initialCondition(rho,v);
end