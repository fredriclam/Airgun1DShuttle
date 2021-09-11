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

function plotPortAndMass(fullState, soln, metadata)

% Pull model data
t = [fullState.t];
pS = [fullState.portStates];
sS = [fullState.shuttleStates];
shuttlePosition = [sS.shuttle_position];
massRate = [pS.massFlowPort];

[~, funcs] = agtools.sampleSignature( ...
    soln, metadata);

global_format = @() set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);

nexttile(1, [2,1])
plot(1e3*t, shuttlePosition, '-', ...
     'LineWidth', 2)
xlim([0,40])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel ('$\xi$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])
yL = ylim();
hold on

nexttile(2, [2,3])
plot(1e3*t, shuttlePosition, '-', ...
     'LineWidth', 2)
xlim([0,300])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
ylim(yL);
hold on
% legend({'Bubble inlet pressure', 'Bubble average pressure'}, ...
%         'Interpreter', 'latex')

nexttile(9, [2,1])
plot(1e3*t, massRate, '-', ...
     'LineWidth', 2)
xlim([0,40])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel ('$\dot{m}$ [kg/s]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])
yL = ylim();
hold on

nexttile(10, [2,3])
plot(1e3*t, massRate, '-', ...
     'LineWidth', 2)
xlim([0,300])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
ylim(yL);
hold on

nexttile(17, [3,1])
plot(1e3*t(1:5:15000), funcs.VDotDotFn(t(1:5:15000)), '-', ...
     'LineWidth', 2);
xlim([0,40])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel ('$\ddot{V}$ [$\mathrm{m}^3$/s]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',4*[0.01, 0.025])
ylim([-2000, 6000])
yL = ylim();
hold on

nexttile(18, [3,3])
plot(1e3*t(1:100:end), funcs.VDotDotFn(t(1:100:end)), '-', ...
     'LineWidth', 2);
xlim([0,300])
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
ylim(yL);
hold on