% Test bubble model with decoupled reservoir
% The purpose of this script is to try out some modifications to add
% turbulence or non-equilibrium effects in the bubble to match the pressure
% sensor data: why is the measured pressure so much lower than the
% thermal equilibrium pressure?
% 
% Needs metadata from a run of airgunShuttleDeploy.

% Get physConst from running the main model
physConst = metadata.discretization.physConst;

V0 = 600*0.0000163871; % [cui -> m^3]
R0 = (3*V0/(4*pi))^(1/3);
Rdot0 = 0;
p_a = 1000*6894.76;
T_a = 300;
rho_a = p_a / (physConst.Q * T_a);
e_a = rho_a * physConst.c_v * T_a;
v_a = sqrt(physConst.gamma * physConst.Q * T_a);

p_inf = 1e5 + 9.8*1000*10; % hydrostatic 10m

rho0 = p_inf / (physConst.Q * T_a);
e0 = rho0 * physConst.c_v * T_a;
m0 = V0 * rho0;
E0 = V0 * e0;
K0 = 0;

bubble0 = [R0; Rdot0; m0; E0; K0];

testRHS = @(t,y,tuningFactor) bubbleRHSTest( ...
    t, y, rho_a, v_a, e_a, p_a, ...
    metadata.discretization.physConst.APortTotal .* (t < 0.1 & t >= 0), ...
    physConst, 'single', tuningFactor);

%% Run ode: uniform pressure p ~ (r/R)^a, a = 1
figure(9); clf;
xLims = [0, 2];
KErelaxationrate = 1/0.0001;
odesoln = ode45(@(t,y)testRHS(t,y,KErelaxationrate), [0, 2], bubble0);
t = linspace(0,2,100000);
[y, dy] = deval(odesoln, t);
% Compute ddot(V)
Vddot = 4*pi*(...
    2*y(1,:) .* y(2,:).^2 ...
    + y(1,:).^2 .* dy(2,:));
% Plot R
subplot(4,1,1);
plot(t,y(1,:));
xlim(xLims)
ylabel 'R [m]'
xlabel 't [s]'
hold on
% Compute rise time
[val, ind] = max(Vddot); 
risetime = t(ind);
% Compute acoustic timescale
subplot(4,1,2);
c = sqrt(y(4,:) ./ y(3,:) * physConst.Q / physConst.c_v * physConst.gamma);
tScale = y(1,:)./c;
plot(t,tScale);
ylabel 'Timescale [s]'
hold on
plot([t(1) t(end)], risetime*[1, 1], '--');
legend({'Acoustic cross-bubble timescale', 'Rise time'})
xlim(xLims)
hold on
% Plot Vddot
subplot(4,1,3);
plot(t,Vddot);
xlim([0, 0.1])
hold on
xlabel 't [s]'; ylabel ('$\ddot{V}$ [m/{s${}^2$}]','Interpreter','latex')
% Compute pressure
V = 4*pi/3*y(1,:).^3;
p = (physConst.gamma-1) * y(4,:) ./ V;
% Plot pressure
subplot(4,1,4);
plot(t,p);
hold on
plot([t(1) t(end)], physConst.p_inf*[1, 1], '--');
xlim([0, 0.1])
ylabel 'p'

% Tuning factor a = -1
hold on
KErelaxationrate = 1/(50.00);
odesoln = ode45(@(t,y)testRHS(t,y,KErelaxationrate), [0, 2], bubble0);
[y, dy] = deval(odesoln, t);
% Compute ddot(V)
Vddot = 4*pi*(...
    2*y(1,:) .* y(2,:).^2 ...
    + y(1,:).^2 .* dy(2,:));
subplot(4,1,1);
plot(t,y(1,:),'m');
subplot(4,1,3);
plot(t,Vddot, 'm');
xlabel 't [s]'; ylabel ('$\ddot{V}$ [m/{s${}^2$}]','Interpreter','latex')
subplot(4,1,4);
V = 4*pi/3*y(1,:).^3;
p = (physConst.gamma-1) * y(4,:) ./ V;
semilogy(t,p,'m');
xlim([0, 0.1])
% ylim([1e5, 1e8])
ylabel 'p'

% % nu = 50*1e-6
% hold on
% odesoln = ode45(@(t,y)testRHS(t,y,5000000), [0, 5], bubble0);
% t = linspace(0,5,1000);
% [y, dy] = deval(odesoln, t);
% % Compute ddot(V)
% Vddot = 4*pi*(...
%     2*y(1,:) .* y(2,:).^2 ...
%     + y(1,:).^2 .* dy(2,:));
% plot(t,Vddot, 'k-.');
% xlabel 't [s]'; ylabel ('$\ddot{V}$ [m/{s${}^2$}]','Interpreter','latex')

function [dy, dQdt, workrate, dEin] = bubbleRHSTest( ...
    t, y, rho_a, v_a, e_a, p_a, A, physConst, bubbleModel, KErelaxationrate)
    R    = y(1);
    Rdot = y(2);
    m    = y(3);
    E    = y(4);
    K = y(5);

    Q       = physConst.Q;
    c_v     = physConst.c_v;
    p_inf   = physConst.p_inf;
    rho_inf = physConst.rho_inf;
    gama    = physConst.gamma;
    c_inf   = physConst.c_inf;

    if strcmpi('quad', bubbleModel)    
        % Surface area and volume factors for hemisphere
        hemisphereFactor = 0.5;
        % Quad bubble rate factor
        rateFactor = 1/4;
    elseif strcmpi('single', bubbleModel)
        hemisphereFactor = 1.0;
        rateFactor = 1;
    end
    
    V = hemisphereFactor * (4/3*pi*R^3);
    Vdot = hemisphereFactor * (4*pi*R^2*Rdot);
    p = E*(gama-1)/V;

    kappa=4000;
    M = 10;
    %M = 25;
    T_inf = 273;
    cv=718;
    Tb = E/(cv*m);
    dQdt = hemisphereFactor*4*pi*R^2*M*kappa*(Tb-T_inf);

    dR = Rdot;
    %dE = A*(e_a + p_a)*v_a - p*Vdot;
    %dE = A*(e_a + p_a)*v_a - p*Vdot - dQdt;
    
    % add turbulent mechanical energy dissipation
    C = 0;
    deltaP = C*rho_inf*abs(Rdot)*Rdot;
    dE = rateFactor*A*(e_a + p_a)*v_a - p*Vdot - dQdt ...
         - hemisphereFactor*4*pi*R^2*Rdot*deltaP;
    workrate = p*Vdot;
    dEin = rateFactor*A*(e_a + p_a)*v_a;
    
    % NEW MECHANICS, replacing above lines
%     KErelaxationrate = 1/(0.0001);
    % Partitioned ratios in choked flow
    kpartition = 64 * gama/2;
    epartition = 1/(gama-1);
    partsize = (kpartition + epartition);
    kpartition = kpartition / partsize;
    epartition = epartition / partsize;
    dK = kpartition*dEin - KErelaxationrate*K;
    dE = epartition*dEin + KErelaxationrate*K - p*Vdot - dQdt ...
         - hemisphereFactor*4*pi*R^2*Rdot*deltaP;

    dpdt = (gama-1)*(dE*V-Vdot*E)/V^2;

    %dRdot = 1/R*((p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt - 3/2*Rdot^2);
    %b = 10; 
    %alpha = b*abs(Rdot); %abs(Rdot);
    b = 0;
    alpha=0*0.8; %b*abs(Rdot); %10;
    nu = 0*1e-6;
    % Nonequilibrium factor assuming power-law pressure solution
    a = 0;
    rarefactionFactor = (a + 3) / 3;
    dRdot = 1/R*((rarefactionFactor*p-p_inf)/rho_inf + R/(rho_inf*c_inf)*dpdt ...
        - 3/2*Rdot^2 - alpha*Rdot ...
        - 4 * nu * Rdot / R); % correction from Langhammer and Landro (1996)
    
    dm = rateFactor*A*rho_a*v_a - 0 * 4*pi*R^2 * Rdot * (m/V);

    dy = [dR; dRdot; dm; dE; dK];
end