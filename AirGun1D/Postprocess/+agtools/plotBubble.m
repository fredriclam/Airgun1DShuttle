%% Plot firing chamber wall state as a function of t
% Plots the specified state at the closed end
% Inputs:
%   state: Full state struct-vector (from compute_full_state)
%   varname: name of field (p, u, T, ...) to plot
%   (optional) rawData: struct(t:time axis, y:data axis) vector with
%   field data objects
% Outputs:
%   t: Extracted t-axis for model
%   field_L: Extracted field at wall for model
% Additional:
%   Plots to current figure.

function plotBubble(fullState, soln, metadata)

% Pull model data
t = [fullState.t];
pS = [fullState.portStates];
bS = [fullState.bubbleStates];
massRate = [pS.massFlowPort];
gamma = metadata.discretization.physConst.gamma;
bubbleMass = [bS.m];

[~, funcs] = agtools.sampleSignature( ...
    soln, metadata);

% Get upstream pressure
pPort = [pS.pPort];
% Compute pressure at inlet of bubble
pExpansionRatio =  ...
    pressureMachFunction(gamma,1) ./ ...
    pressureMachFunction(gamma, [pS.MPort]);
% Pressure if choked
pBubblePortChoked = pPort .* pExpansionRatio;
% Compute thermodynamic pressure
pBubbleThermo = [bS.p];

% Get inlet pressure
pBubbleInlet = pBubbleThermo;
pBubbleInlet(pBubblePortChoked > pBubbleThermo) = ...
    pBubblePortChoked(pBubblePortChoked > pBubbleThermo);

global_format = @() set(gca, 'FontSize', 14, ...
    ...'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);

tL = tiledlayout(3,3);

nexttile(tL, 1, [1,1])
plotStep = 1;
% plot(tAxis(1:plotStep:end), pPort(1:plotStep:end)/1e6, 'k')
plot(1e3*t(1:plotStep:5000), pBubbleInlet(1:plotStep:5000)/1e6, '--', ...
     'LineWidth', 1, 'Color', [123, 31, 21]/255)
hold on
plot(1e3*t(1:plotStep:5000), pBubbleThermo(1:plotStep:5000)/1e6, ...
     'k', 'LineWidth', 1)
hold off
xlim([0,10])
xlabel('{\it{t}} (ms)', 'FontSize', 14)
ylabel ('{\it{p}}_b (MPa)', 'FontSize', 14)
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel ('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])
yL = ylim();

nexttile(tL, 2, [1,2])
plot(1e3*t(1:100:end), pBubbleInlet(1:100:end)/1e6, '--', ...
     'LineWidth', 1, 'Color', [123, 31, 21]/255)
hold on
plot(1e3*t(1:100:end), pBubbleThermo(1:100:end)/1e6, ...
     'k', 'LineWidth', 1)
hold off
xlim([0,300])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('{\it{t}} (ms)', 'FontSize', 14)
% ylabel ('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
ylim(yL);
legend({'Bubble inlet pressure', 'Bubble average pressure'}, ...
        'Interpreter', 'latex')

nexttile(tL, 4, [1,1])
plot(1e3*t(1:plotStep:5000), bubbleMass(1:plotStep:5000), 'k-', ...
     'LineWidth', 1)
xlim([0,10])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('{\it{t}} (ms)', 'FontSize', 14)
% ylabel ('$m$ [kg]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel ('{\it{m}} (kg)', 'FontSize', 14)

set(gca,'XMinorTick','on', 'YMinorTick', 'on', 'LineWidth', 1)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])

nexttile(tL, 5, [1,2])
plot(1e3*t(1:100:end), bubbleMass(1:100:end), 'k-', ...
     'LineWidth', 1)
xlim([0,300])
ylim([0, 30])
xlabel('{\it{t}} (ms)', 'FontSize', 14)
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel ('$m$ [kg]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()

 
nexttile(tL, 7, [1,1])
plot(1e3*t(1:plotStep:5000), funcs.VDotDotFn(t(1:plotStep:5000)), 'k-', ...
     'LineWidth', 1)
xlim([0,10])
xlabel('{\it{t}} (ms)', 'FontSize', 14)
ylabel ('$\ddot{V}$ (m${}^3$/s${}^2$)', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])

nexttile(tL, 8, [1,2])
plot(1e3*t(1:100:end), funcs.VDotDotFn(t(1:100:end)), 'k-', ...
     'LineWidth', 1)
xlim([0,300])
xlabel('{\it{t}} (ms)', 'FontSize', 14)
ylabel ('$\ddot{V}$ (m${}^3$/s${}^2$)', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
ylim([-3e3, 6e3])
set(gca,'YTick',[-3000, 0, 3000, 6000])

% plot(1e3*t, massRate, 'k', 'LineWidth', 1);
% 
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex', 'FontSize', 14)
% set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
%     'XMinorTick', 'on', 'YMinorTick', 'on');