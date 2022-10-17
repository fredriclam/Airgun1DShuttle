function tL = plotBubble_light(fullState, soln, metadata, tL, lineColor)
if nargin < 4
    tL = tiledlayout(2,1);
end
if nargin < 5
    lineColor = [0 0 0];
end

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

plotStep = 10;
global_format = @() set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);

nexttile(tL, 1, [1,1])
plot(1e3*t(1:plotStep:7000), bubbleMass(1:plotStep:7000), 'k-', ...
     'LineWidth', 1, 'Color', lineColor)
xlim([0,15])
xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 14)
ylabel ('{\it{m}} (kg)', 'Interpreter', 'tex', 'FontSize', 14)
set(gca,'XMinorTick','on', 'YMinorTick', 'on', 'LineWidth', 1)
global_format()
set(gca,'TickLength',1.5*[0.01, 0.025])
hold on

nexttile(tL, 2, [1,1])
plot(1e3*t(1:plotStep:7000), funcs.VDotDotFn(t(1:plotStep:7000)), 'k-', ...
     'LineWidth', 1, 'Color', lineColor)
xlim([0,15])
xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 14)
ylabel ('$\ddot{V}$ (m${}^3$/s${}^2$)', 'Interpreter', 'latex', 'FontSize', 14)
global_format()
set(gca,'TickLength',1.5*[0.01, 0.025])
hold on
