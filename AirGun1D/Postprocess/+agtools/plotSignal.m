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