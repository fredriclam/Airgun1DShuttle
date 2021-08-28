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

function [tSample, p, dists, funcs] = plotSignal(...
    fullState, soln, metadata, hydrophoneDepth, ...
    lateralSeparation, tSample, color, linestyle)

% Set default sampling vector if not provided
if nargin < 6
    tSample = linspace(0,5,7500);
end

% Set default linecolor
if nargin < 7
    color = [0, 0, 0];
end

% Set defualt linestyle
if nargin < 8
    linestyle = '--';
end

% Pull model data
t = [fullState.t];
pS = [fullState.portStates];
bS = [fullState.bubbleStates];
massRate = [pS.massFlowPort];
gamma = metadata.discretization.physConst.gamma;
bubbleMass = [bS.m];

[p, funcs] = agtools.sampleSignature( ...
    soln, metadata, tSample, hydrophoneDepth, lateralSeparation);

% Estimate of effective floor depth for reflection
floorDepth = 27;
depth = metadata.paramAirgun.airgunDepth;
% Compute distances from source, surface ghost, floor ghost
r1 = norm([lateralSeparation, depth-hydrophoneDepth]);
r2 = norm([lateralSeparation, depth+hydrophoneDepth]);%
r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepth]);
dists = struct('r1', r1, 'r2', r2, 'r3', r3);

global_format = @() set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);

% tL = tiledlayout(3,4);
% nexttile(tL, 1, [1,1])
plotStep = 1;
plot(1e3*tSample, p*r1/1e5, linestyle, ...
     'LineWidth', 1, 'Color', color)
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel ('$\Delta p \cdot r_1$ [bar $\cdot$ m]', ...
        'Interpreter', 'latex', 'FontSize', 14)
global_format();

% plot(1e3*t(1:plotStep:5000), pBubbleThermo(1:plotStep:5000)/1e6, ...
%      'k', 'LineWidth', 1)
% hold off
% xlim([0,10])
% global_format()
% set(gca,'TickLength',4*[0.01, 0.025])
% yL = ylim();
% 
% nexttile(tL, 2, [1,3])
% plot(1e3*t(1:100:end), pBubbleInlet(1:100:end)/1e6, '--', ...
%      'LineWidth', 1, 'Color', [123, 31, 21]/255)
% hold on
% plot(1e3*t(1:100:end), pBubbleThermo(1:100:end)/1e6, ...
%      'k', 'LineWidth', 1)
% hold off
% xlim([0,300])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% % ylabel ('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
% global_format()
% ylim(yL);
% legend({'Bubble inlet pressure', 'Bubble average pressure'}, ...
%         'Interpreter', 'latex')
% 
% nexttile(tL, 5, [1,1])
% plot(1e3*t(1:plotStep:5000), bubbleMass(1:plotStep:5000), 'k-', ...
%      'LineWidth', 1)
% xlim([0,10])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel ('$m$ [kg]', 'Interpreter', 'latex', 'FontSize', 14)
% set(gca,'XMinorTick','on', 'YMinorTick', 'on', 'LineWidth', 1)
% global_format()
% set(gca,'TickLength',4*[0.01, 0.025])
% 
% nexttile(tL, 6, [1,3])
% plot(1e3*t(1:100:end), bubbleMass(1:100:end), 'k-', ...
%      'LineWidth', 1)
% xlim([0,300])
% ylim([0, 30])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% % ylabel ('$m$ [kg]', 'Interpreter', 'latex', 'FontSize', 14)
% global_format()
% 
%  
% nexttile(tL, 9, [1,1])
% plot(1e3*t(1:plotStep:5000), funcs.VDotDotFn(t(1:plotStep:5000)), 'k-', ...
%      'LineWidth', 1)
% xlim([0,10])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel ('$\ddot{V}$ [m${}^3$/s${}^2$]', 'Interpreter', 'latex', 'FontSize', 14)
% global_format()
% set(gca,'TickLength',4*[0.01, 0.025])
% 
% nexttile(tL, 10, [1,3])
% plot(1e3*t(1:100:end), funcs.VDotDotFn(t(1:100:end)), 'k-', ...
%      'LineWidth', 1)
% xlim([0,300])
% xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel ('$\ddot{V}$ [m${}^3$/s${}^2$]', 'Interpreter', 'latex', 'FontSize', 14)
% global_format()
% ylim([-2e3, 6e3])
