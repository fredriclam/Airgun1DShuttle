% Generates Figure 4 (timescales)

%% Simulate input
solution = solution_reference;
metadata = metadata_reference;
hydrophoneDepth = 9;
lateralSeparation = 6;

%% Unpack system states
bS = [fullState.bubbleStates];
eDS = [fullState.eulerDomainStates];
pS = [fullState.portStates];
sS = [fullState.shuttleStates];

%% Panel (a)
figure(4);
tL = tiledlayout(9,1);
nexttile(tL, 1, [3,1]);

semilogx(solution.soln.x, [bS.p]/1e6, 'LineWidth', 1);
hold on
semilogx(solution.soln.x, [pS.pPort]/1e6, 'LineWidth', 1);
extentx = [solution.soln.x(2) solution.soln.x(end)];
semilogx(extentx, 1000*0.00689476*ones(size(extentx)), 'k--', ...
    'LineWidth', 0.5);
hold off
xlim([1e-4, 1])
ylim([0, 7])

set(gca, 'FontSize', 12, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLabelInterpreter', 'tex')
xlabel("{\it{t}} (s)")
ylabel("{\it{p}} (MPa)")
grid on

legend(["Bubble pressure", "{\it{p}} ({\it{x}} = 0)"])

%% Panel (b)
nexttile(tL, 4, [2,1]);

semilogx(solution.soln.x, [sS.shuttle_position], 'k', 'LineWidth', 1);
xlim([1e-4, 1])

set(gca, 'FontSize', 12, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLabelInterpreter', 'tex')
xlabel("{\it{t}} (s)")
ylabel("\xi (m)")
grid on

%% Panel (c)
nexttile(tL, 6, [2,1]);

tBubble = solution.bubbleContinuationTime;
RBubble = solution.bubbleContinuationState(1,:);
loglog(tBubble, 4/3*pi*RBubble.^3, 'k', 'LineWidth', 1);
xlim([1e-4, 1])
ylim([1e-5, 1e2])

set(gca, 'FontSize', 12, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLabelInterpreter', 'tex')
xlabel("{\it{t}} (s)")
ylabel("{\it{V}} (m^3)")
grid on

%% Panel (d)
nexttile(tL, 8, [2,1]);

c_inf = metadata.discretization.physConst.c_inf;

tSample = linspace(0,1,10000);
% Gobble plot output
figure(994);
disp("Computing acoustic pressure.")
[tSample, p, dists, funcs] = agtools.plotSignal(...
    fullState, solution, metadata, hydrophoneDepth, ...
    lateralSeparation, tSample);%, color, linestyle)
disp("Done acoustic pressure.")
drawnow

figure(4);
drawnow
semilogx(tSample - dists.r1 / c_inf, p/1e3, 'k', 'LineWidth', 1);
xlim([1e-4, 1])
ylim([-40, 80])

set(gca, 'FontSize', 12, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLabelInterpreter', 'tex', ...
    'YTick', [-40, 20, 80])
xlabel("{\it{t}} - {\it{r}}_1 / {\it{c}}_\infty (s)")
ylabel("\Delta {\it{p}} (kPa)")
grid on

%% Set figure position
set(gcf, 'position', [680   250   566   727])

%% Tool: port area history with superimposed cross-sectional area
figure(104); clf
semilogx(solution.soln.x, [pS.APortExposed])
hold on
semilogx(solution.soln.x, ...
    metadata.paramAirgun.airgunCrossSecAreaSqInch*0.0254^2*ones(size([pS.APortExposed])))
xlim([1e-4, 1])