% Test bubble model with decoupled reservoir
% The purpose of this script is to try out some modifications to add
% turbulence or non-equilibrium effects in the bubble to match the pressure
% sensor data: why is the measured pressure so much lower than the
% thermal equilibrium pressure?
% 
% Needs metadata from a run of airgunShuttleDeploy.

% Get physConst from running the main model
physConst = metadata_reference.discretization.physConst;

V0 = 1*0.0000163871; % [cui -> m^3]
R0 = (3*V0/(4*pi))^(1/3);
Rdot0 = 0;
p_a = 1000*6894.76;
T_a = 288;
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
%% Extract solved coupled-model data
tAxis = [fullState.t];
pS = [fullState.portStates];
massrate = [pS.massFlowPort];

rhoPort = [pS.rhoPort];
vPort = [pS.velocityPort];
eTotalPort = [pS.eTotalPort];
pPort = [pS.pPort];

funcify = @(t, data) pchip(tAxis, data, t);

KErelaxationrate = 1/0.010;
testRHS = @(t,y) ke_bubble( ...
    t, y, ...
    funcify(t, rhoPort), ...
    funcify(t, vPort), ...
    funcify(t, eTotalPort), ...
    funcify(t, pPort), ...
    0.8 * physConst.APortTotal ...
       .* (t < 0.300 & t >= 0), ...
    physConst, 'single', KErelaxationrate);

%% Run ode: uniform pressure p ~ (r/R)^a, a = 1
odesoln = ode45(@(t,y) testRHS(t,y), [0, 0.300], bubble0);
disp('ODE done')

%% Plot
figure(999); clf;
tiledlayout('flow')
xLims = [0, 0.3];
t = linspace(0,0.3,10000);

% Compute ddot(V)
[y, dy] = deval(odesoln, t);
Vddot = 4*pi*(...
    2*y(1,:) .* y(2,:).^2 ...
    + y(1,:).^2 .* dy(2,:));

% Get ref solution
[yReference, dyReference] = deval( ...
    soltest_.soln, t);
% Slice for bubble states only
yReference = yReference(end-10:end-6,:);
dyReference = dyReference(end-10:end-6,:);

tReference = t; % solution_reference.soln.x;
% Compute ddot(V)
VddotReference = 4*pi*(...
    2*yReference(1,:) .* yReference(2,:).^2 ...
    + yReference(1,:).^2 .* dyReference(2,:));

% Plot R
nexttile;
plot(t,y(1,:));
hold on
plot(tReference,yReference(1,:));
hold off
xlim(xLims)
ylabel 'R [m]'
xlabel 't [s]'
legend({'k-e', 'ref'})

if false
    % Compute rise time
    [val, ind] = max(Vddot); %#ok<UNRCH>
    risetime = t(ind);
    % Compute acoustic timescale
    nexttile;
    c = sqrt(y(4,:) ./ y(3,:) * physConst.Q / physConst.c_v * physConst.gamma);
    tScale = y(1,:)./c;
    plot(t,tScale);
    ylabel 'Timescale [s]'
    hold on
    plot([t(1) t(end)], risetime*[1, 1], '--');
    legend({'Acoustic cross-bubble timescale', 'Rise time'})
    xlim(xLims)
    hold on
end

% Plot Vddot
nexttile;
plot(t,Vddot);
xlim(xLims)
xlabel 't [s]'; ylabel ('$\ddot{V}$ [m/{s${}^2$}]','Interpreter','latex')
% Compute pressure
V = 4*pi/3*y(1,:).^3;
p = (physConst.gamma-1) * y(4,:) ./ V;

VReference = 4*pi/3*yReference(1,:).^3;
pReference = (physConst.gamma-1) * yReference(4,:) ./ VReference;

hold on
plot(tReference,VddotReference);
legend({'k-e', 'ref'})

% Plot pressure
nexttile;
plot(t,p);
hold on
plot(tReference,pReference);
hold off

% plot([t(1) t(end)], physConst.p_inf*[1, 1], '--');
xlim([0, 0.1])
ylabel 'p'
legend({'k-e', 'ref'})
xlim([0, 1e-2])
