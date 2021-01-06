function [icAirgun, icBubble, t0] = icAirgunShocktube(t0 , x0, A, physConst)
    default_arg('t0',[]);
    % A -- area of airgun port

    gamma = physConst.gamma;  %
    R_G   = physConst.R_G;    % Specific gas constant
    Tinf  = physConst.Tinf;   % Temperature assumed constant throughout the system
    c_v   = physConst.c_v;    %


    n = 7;
    B = 303.975e6;   %??????????????????????

    %% Region 1  (outside)
    p1 = physConst.p_inf; %Pa
    v1 = 0;
    rho1 = physConst.rho_inf; %kg/m^3
    c1 = sqrt(n*(p1+B)/rho1); % 1459 m/s;

    %% Region 4  (inside)
    p4 = physConst.p0a; %Pa
    v4 = 0;
    c4 = sqrt(gamma*R_G*Tinf);  %% 335m/s
    rho4 = p4/(R_G*Tinf);  % 124 kg/m^3


    %% Region 2  (outside)
    rho2 = @(p2) rho1*((p2+B)/(p1+B))^(1/n);
    v2   = @(p2) sqrt((p2-p1)*(rho2(p2)-rho1)/(rho1*rho2(p2)));
    % p2 is the solution to the equation:
    f = @(p2) p4 - p2*(1 + (gamma - 1)/(2*c4)*(v4-v2(p2)))^(2*gamma/(1-gamma));
    p2 = fzero(f, p1);
    rho2 = rho2(p2);
    v2 = v2(p2);

    %% Region 3  (inside)
    v3 = v2;
    p3 = p2;
    c3 = (gamma-1)/2 * (v4 - v3 + 2*c4/(gamma - 1));
    rho3 = gamma*p3/c3^2;


    if isempty(t0)
        t0 = sqrt(A)/v2; % Choose starting time so that bubble is 'square'??
    end


    %% Expansion fan
    v_exp_upstream  = v4 - c4;
    v_exp_downstream = v3 - c3;

    v_e   = @(x) 2/(gamma + 1) * ((x-x0)/t0 + (gamma -1)/2*v4 + c4);
    c_e   = @(x) v_e(x)-(x-x0)/t0;
    p_e   = @(x) p4*(c_e(x)/c4).^(2*gamma/(gamma - 1));
    rho_e = @(x) gamma*p_e(x)./c_e(x).^2;


    %% Airgun initial conditions
    x_up   = x0 + v_exp_upstream*t0;
    x_down = x0 + v_exp_downstream*t0;

    if abs(v_exp_upstream*t0) > physConst.L;
        warning('The expansion fan traveled longer that the length of the airgun before the first timestep.\nabs(v_exp_upstream*t0)=%f',abs(v_exp_upstream*t0))
    end

    v   = @(x)(x<=x_up)*v4   + (x>x_up & x<=x_down).*v_e(x)    +  (x>x_down)*v3  ;
    rho = @(x)(x<=x_up)*rho4 + (x>x_up & x<=x_down).*rho_e(x)  +  (x>x_down)*rho3;
    p   = @(x)(x<=x_up)*p4   + (x>x_up & x<=x_down).*p_e(x)    +  (x>x_down)*p3  ;
    e   = @(x)p(x)/(gamma-1)+1/2*rho(x).*v(x).^2; %???????????????????????

    icAirgun.rho0 = rho;
    icAirgun.rv0  = @(x)rho(x).*v(x);
    icAirgun.e0   = e;
    icAirgun.p0   = p;


    %% Bubble initial conditions
    V = A*v2*t0;
    p = p2;

    icBubble.R = (3/(4*pi) * V)^(1/3);
    icBubble.Rdot = v2;
    icBubble.m = rho2*V;
    icBubble.E = V * p/(gamma - 1); %????????????????????????????;

end