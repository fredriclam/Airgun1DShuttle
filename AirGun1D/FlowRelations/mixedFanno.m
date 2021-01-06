% Area change plus Fanno flow integration
%
% Taking perfect gas, mass cons., streamwise momentum, and energy
% conservation with averaged friction coefficient C_f and arbitrary area
% change function (and no stagnation energy change)

clear;
clc; figure(5); clf;

%% Constants and parameters
gamma = 1.4;
Q = 287;
cp = Q * gamma / (gamma - 1);
% cp = 1000;

LDRatio = 10;  % Effective length divided by average hydraulic diameter
A1 = 0.0077;
A2 = A1/1;

A = @(x) x * (A2-A1) + A1; % Linear interpolated area function, nondim x
dAdx = @(x) 1 * (A2-A1);   % Area profile derivative

%% Initial
T_initial = 300;
p_initial = 800 * 6895;
rho_initial = p_initial / (Q * T_initial);
u2_initial = (10)^2;
M_initial = sqrt(u2_initial)/sqrt(gamma * Q * T_initial);
% Initial stagnation values
T0 = T_initial + 0.5 * u2_initial / cp;
p0 = p_initial / pressureMachFunction(gamma, M_initial);
% Initial state vector
y0 = [rho_initial; T_initial; p_initial; u2_initial];

%% Analytical zero-friction limit
machAreaFunction = precomputeMachAreaFunction(gamma);
xVector = linspace(0,1,1000);
machVector = machAreaFunction(A(xVector) / A1 * areaMachFunction(gamma, M_initial));
temperatureVector = temperatureMachFunction(gamma, machVector) * T0;
pressureVector = pressureMachFunction(gamma, machVector) * p0;
rhoVector = pressureVector ./ (Q * temperatureVector);

%% Base layer plotting
subplot(2,2,1);
plot(xVector, rhoVector,'-k','LineWidth',1.0); hold on;
ylabel('$\rho$ [kg/m${}^3$]', 'Interpreter', 'latex')

subplot(2,2,2);
plot(xVector, temperatureVector,'-k','LineWidth',1.0); hold on;
ylabel('$T$ [K]', 'Interpreter', 'latex')

subplot(2,2,3);
plot(xVector, pressureVector,'-k','LineWidth',1.0); hold on;
ylabel('$p$ [Pa]', 'Interpreter', 'latex')

subplot(2,2,4);
plot(xVector, sqrt(2 * (T0 - temperatureVector) * cp),'-k','LineWidth',1.0); hold on;
ylabel('$u$ [m/s]', 'Interpreter', 'latex')

for i = 1:4
    subplot(2,2,i);
    xlabel('$x = x_\mathrm{dim}/L$', 'Interpreter', 'latex')
end

disp('Analytic done.');

%% Friction addition
overPlot(0.0, LDRatio, A, dAdx, gamma, Q, y0);
overPlot(0.1, LDRatio, A, dAdx, gamma, Q, y0); % Upper bound on incompressible Re 1e8 flow
overPlot(1.0, LDRatio, A, dAdx, gamma, Q, y0);
overPlot(5.0, LDRatio, A, dAdx, gamma, Q, y0);
% overPlot(20.0, LDRatio, A, dAdx, gamma, Q, y0);
% overPlot(1000.0, D, A, dAdx, gamma, Q, y0);

subplot(2,2,2);
legend({'$C_f = 0; L/D = 10$ analytic', ...
        '$C_f = 0; L/D = 10$', ...
        '$C_f = 0.1; L/D = 10$', ...
        '$C_f = 1.0; L/D = 10$', ...
        '$C_f = 5.0; L/D = 10$'},'Interpreter', 'latex');
disp('Numerics done.');

%% Definition of flow RHS
function dy = flowEvolution(y, x, Cf, LDRatio, A, dAdx, gamma, Q)
    % Unpack state
    rho = y(1);
    T = y(2);
    p = y(3);
    u2 = y(4); % u^2 state variable

    % Compute helpers
    M = sqrt(u2)/sqrt(gamma * Q * T);
    c1 = 2./ (1 - M.^2);
    f1 = Cf * LDRatio;
    
    % Compute differentials
    drhodx = ...
        c1 * rho * ...
        (0.5 * M.^2 * dAdx(x) / A(x) - ...
        0.25 * gamma * M.^2 * f1);
    dTdx = ...
        c1 * T * ...
        (0.5 * (gamma-1) * M.^2 * dAdx(x) ./ A(x) - ...
        0.25 * gamma * (gamma-1) * M.^4 * f1);
    dpdx = ...
        c1 * p * ...
        (0.5 * gamma * M.^2 * dAdx(x) / A(x) - ...
        gamma * 0.25 * M.^2 * (1 + (gamma-1)*M.^2) * f1);
    du2dx = ...
        c1 * u2 * (-dAdx(x) / A(x) + 0.5 * gamma * M.^2 * f1);
    
    % Pack differential
    dy = [drhodx; dTdx; dpdx; du2dx];
end

% Overlay a plot with friction factor
function overPlot(Cf, LDRatio, A, dAdx, gamma, Q, y0)
    %% Run ODE -- mixed CD-Fanno flow
    soln = ode45(@(x,y) flowEvolution(y, x, Cf, LDRatio, A, dAdx, ...
        gamma, Q), ...
        [0, 1], y0);

    %% Add to plot
    for i = 1:3
        subplot(2,2,i);
        plot(soln.x, soln.y(i,:),'--', 'LineWidth', 1.5);
    end
    subplot(2,2,4);
    plot(soln.x, sqrt(soln.y(4,:)),'--', 'LineWidth', 1.5);
end