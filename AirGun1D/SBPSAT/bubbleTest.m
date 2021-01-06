t0 = 0.001;
A = 0.001;
xlim = {0, 1};
T = 100*t0;

physConst.R_G = 287; %?????????
physConst.c_v = 718; %????? For dry air, is this what goes into the IC? another one for water?
physConst.p_inf = 1e5; % Is this p1 from IC? should it be used for p1 in IC?
physConst.p0a = 1e7;
physConst.rho_inf = 1e3; %????? For water? should this also be used for rho1 in IC?
physConst.Tinf = 280; %K         Temperature assumed constant throughout the system
physConst.gamma = 1.4;
physConst.L = 1;

[icAirgun, icBubble] = icAirgunShocktube(t0, xlim{2}, A, physConst);

y0 = [
    icBubble.R;
    icBubble.Rdot;
    icBubble.m;
    icBubble.E;
];

rho_a = icAirgun.rho0(xlim{2});
v_a = icAirgun.rv0(xlim{2})/rho_a;
e_a = icAirgun.e0(xlim{2});
p_a = icAirgun.p0(xlim{2});

f = @(t,y)bubbleRHS(y, rho_a, v_a, e_a, p_a, A, physConst);
% f = @(t,y)bubbleRHS(y, rho_a, 0, e_a, p_a, A, physConst);

[t,y] = ode15s(f,[t0 T], y0,odeset('OutputFcn',@odeplot));
size(y)
plot(t,y(:,1));