%% Rayleigh-Plesset asymptotics
epsilon = 1e-7;
rho_a = 50;
v_a = 300;
T_A = 288;
c_v = 718;
gamma = 1.4;
c_p = gamma * c_v;
c_vT = c_v*T_A;
c_pT = c_p*300;
e_a = rho_a * c_vT;
p_a = 5e6;
A = 0.0465;

y = [epsilon,
     epsilon,
     4/3*pi*epsilon^3*rho_a,
     4/3*pi*epsilon^3*e_a];
physConst = struct('Q', 287.06, 'c_v', c_v, 'p_inf', 1e5, 'rho_inf', 1e3 ...
    ,'gamma', gamma, 'c_inf', 1482);
dy = @(t, y) bubbleRHS_Unmodified(y, rho_a, v_a, e_a, p_a, A, physConst);
out = ode45(dy, [0,5], y);

%% Plot figure with epsilon noise
figure(1);
clf
subplot(1,2,1);
loglog(out.x, out.y(1,:))
hold on

% Compute power law prefactor
mdot = rho_a * v_a * A;
rho_inf = physConst.rho_inf;
beta = 3*(gamma - 1)/(4 * pi * rho_inf) * mdot * c_pT;
prefactor = (8 * beta)^(1/5);
prefactor = (beta / 0.3)^(1/5);

plot(logspace(-15,1,100),epsilon+prefactor*logspace(-15,1,100).^(3/5))
plot([1e-20, 1e2],1e-2*[1e-20, 1e2].^(1/3), 'k--')
figure(gcf);
ylim([1e-8, 1e1])

%% Extrapolation on figure
plot(logspace(-15,1,100),prefactor*logspace(-15,1,100).^(3/5))
% plot(logspace(-15,1,100),0.8*prefactor*logspace(-15,1,100).^(3/5))
legend({'Bubble model (mod RP)', '3/5 law w/ eps noise', ...
    '1/3 reference', '3/5 law extrapolation'}, 'Location', 'best')
title ('Rayleigh-Plesset, sans correction')

%% Linear scales
subplot(1,2,2);
timescale = 0.010;
plot(out.x(out.x < timescale), out.y(1, out.x < timescale))
hold on
plot(out.x(out.x < timescale), prefactor * out.x(out.x < timescale).^(3/5))
plot(out.x(out.x < timescale), 0.9 * prefactor * out.x(out.x < timescale).^(3/5))
legend({'Bubble model (mod RP)','3/5 law prediction', '3/5 with fudging'})
title ('Rayleigh-Plesset, sans correction')
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1:2), 1200, 420]);




%% Modified MP comparison
epsilon = 1e-7;
rho_a = 50;
v_a = 300;
c_vT = 700*300;
c_pT = 1000*300;
e_a = rho_a * c_vT;
p_a = 5e6;
A = 0.0465;

y = [epsilon,
     epsilon,
     4/3*pi*epsilon^3*rho_a,
     4/3*pi*epsilon^3*e_a];
physConst = struct('Q', 287, 'c_v', 700, 'p_inf', 1e5, 'rho_inf', 1e3 ...
    ,'gamma', 1.4, 'c_inf', 1e3);
dy = @(t, y) bubbleRHS(y, rho_a, v_a, e_a, p_a, A, physConst);
out = ode45(dy, [0,5], y);

%% Plot figure with epsilon noise
figure(2);
clf
subplot(1,2,1);
loglog(out.x, out.y(1,:))
hold on

% Compute power law prefactor
mdot = rho_a * v_a * A;
rho_inf = 1000;
gamma = 1.4;
beta = 3*(gamma - 1)/(4 * pi * rho_inf) * mdot * c_pT;
prefactor = (8 * beta)^(1/5);

plot(logspace(-15,1,100),epsilon+prefactor*logspace(-15,1,100).^(3/5))
plot([1e-20, 1e2],1e-2*[1e-20, 1e2].^(1/3), 'k--')
figure(gcf);
ylim([1e-8, 1e1])

%% Extrapolation on figure
plot(logspace(-15,1,100),prefactor*logspace(-15,1,100).^(3/5))
plot(logspace(-15,1,100),15*logspace(-15,1,100).^(3/5))
legend({'Bubble model (mod RP)', '3/5 law w/ eps noise', ...
    '1/3 reference', '3/5 law extrapolation', ...
    '3/5 with fudging'}, ...
    'Location', 'best')
title ('Rayleigh-Plesset, with correction')

%% Linear scales
subplot(1,2,2);
timescale = 0.010;
plot(out.x(out.x < timescale), out.y(1, out.x < timescale))
hold on
plot(out.x(out.x < timescale), prefactor * out.x(out.x < timescale).^(3/5))
legend({'Bubble model (mod RP)','3/5 law prediction'})
title ('Rayleigh-Plesset, with correction')
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1:2), 1200, 420]);