function ic = icShuOsher()
    gamma = 1.4;

    wavelength = 1/8;   % To Time   0.178
    x0 = 1/8;

    % wavelength = 2/5;   % To time 0.47
    % x0 = -4;

    rho_l = 3.857143;
    rho_r = @(x)1+0.2*sin(2*pi*1/wavelength*x);

    v_l = 2.629369;
    v_r = 0;

    p_l = 10.3333;
    p_r = 1;


    c = 1; %speed of sound ????????????

    v   = @(x)(x<=x0)*v_l   + (x>x0)*v_r;
    rho = @(x)(x<=x0)*rho_l + (x>x0).*rho_r(x);
    p   = @(x)(x<=x0)*p_l   + (x>x0)*p_r;
    e   = @(x)p(x)/(gamma-1)+1/2*rho(x).*v(x).^2;

    ic.rho0 = rho;
    ic.rv0  = @(x)rho(x).*v(x);
    ic.e0   = e;
    ic.p0   = p;
end

% This test is a hydrodynamic shocktube where the left and the right states are given as follows.
% Left: (ρ=3.857143; Vx= 2.629369, P = 10.33333) Right: (ρ=1 + 0.2 sin(5 π x); Vx=0; P=1).
% Essentially it is a Mach=3 shock interacting with a sine wave in density.


% ρ=3.857143;            u=2.629369;          P=10.3333             when x<1/8
% ρ=1+0.2sin8x;          u=0;                    P=1                      when x≥1/8


% fine grid 192
% coarse 384