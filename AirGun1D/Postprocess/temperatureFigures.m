%% Different airgun temperature figures
% Build figures for T_inf not equal to the temperature at the initial
% condition of the firing chamber (and bubble). The temperature of the
% operating chamber is the same.

% Working script that contains various postprocess routines.
% To be split into multiple files.

%% Header
% Add required dependencies from AirGun1D root.
% The following assumes the user is in the ./AirGun1D folder (otherwise,
% replace the following lines with the relative path to the respective
% folders.
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Global parameters
nx = 40;

%% Bubble comparison
% Background process or parallel (requires O(10 GB) per thread)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(3);
end

futuresBubbleModel(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single', 'extraOptions', ...
      struct('TInitial', 288)));
futuresBubbleModel(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single', 'extraOptions', ...
      struct('TInitial', 288-15)));
futuresBubbleModel(3) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single', 'extraOptions', ...
      struct('TInitial', 288+15)));

%% Run uncoupled model
% futuresComparison = parfeval(pool, @airgunShuttleDeploy, 2, ...
%     nx, false);

%% Reduce data size
% for j = 1:5
%     dataPortAreaRatio(j) = dataReduction(futuresPortAreaRatio(j));
% end
% for j = 1:6
%     dataPressures(j) = dataReduction(futuresPressures(j));
% end

%% Process full state for reference cases

[fullState, caseKeyContext] = ...
        airgunShuttlePostprocess( ...
        futuresBubbleModel(1).OutputArguments{1}, ...
        futuresBubbleModel(1).OutputArguments{2});
% [fullStateQuad, caseKeyContextQuad] = ...
%         airgunShuttlePostprocess( ...
%         futuresBubbleModel(2).OutputArguments{1}, ...
%         futuresBubbleModel(2).OutputArguments{2});
% [fullStateUncoupled] = ...
%         airgunShuttlePostprocess( ...
%         futuresComparison.OutputArguments{1}, ...
%         futuresComparison.OutputArguments{2});


%% debug
% [fullState_, caseKeyContext_] = ...
%         airgunShuttlePostprocess( ...
%         futuresBubbleModel(1).OutputArguments{1}, ...
%         futuresBubbleModel(1).OutputArguments{2});

%% Acoustic pressure signals
% Compose tSample from several segments
tEvent = 0.004104428158096; % t at which first wave arrives
tSample_1 = linspace(0,0.95*tEvent,100);
tSample_2 = linspace(0.95*tEvent,4*tEvent,1000);
tSample_3 = linspace(4*tEvent, 5, 15000);
tSample = [tSample_1(1:end-1), ...
           tSample_2(1:end-1), ...
           tSample_3];

[signals15, tFirst] = sampleSignature( ...
        futuresBubbleModel(1).OutputArguments{1}, ...
        futuresBubbleModel(1).OutputArguments{2}, ...
        tSample);
signalsCold = sampleSignature( ...
        futuresBubbleModel(2).OutputArguments{1}, ...
        futuresBubbleModel(2).OutputArguments{2}, ...
        tSample);
signalsHot = sampleSignature( ...
        futuresBubbleModel(3).OutputArguments{1}, ...
        futuresBubbleModel(3).OutputArguments{2}, ...
        tSample);
    
%% Plot
figure(11); clf;
plot(tSample, signalsCold)
hold on
plot(tSample, signals15)
plot(tSample, signalsHot)
hold off

legend({'0 deg C', '15 deg C', '30 deg C'})

%% Plots using full state
tAxis = [fullState.t];
eDS = [fullState.eulerDomainStates];
pS = [fullState.portStates];
bS = [fullState.bubbleStates];
sS = [fullState.shuttleStates];

% Get upstream pressure
pPort = [pS.pPort];

% Compute pressure at inlet of bubble
pExpansionRatio =  ...
    pressureMachFunction(1.4,1) ./ pressureMachFunction(1.4, [pS.MPort]);
% Pressure if choked
pBubblePortChoked = pPort .* pExpansionRatio;
pBubbleInlet = pBubbleThermo;
pBubbleInlet(pBubblePortChoked > pBubbleThermo) = ...
    pBubblePortChoked(pBubblePortChoked > pBubbleThermo);
% Compute thermodynamic pressure
pBubbleThermo = [bS.p];

figure(36); clf;
global_format = @() ...
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex', ...
        'LineWidth', 0.75);

plotStep = 1;
tL = tiledlayout(9,1);
nexttile(tL, 1, [3,1])
% plot(tAxis(1:plotStep:end), pPort(1:plotStep:end)/1e6, 'k')
plot(tAxis(1:plotStep:end), pBubbleInlet(1:plotStep:end)/1e6, '--', ...
     'LineWidth', 1, 'Color', [123, 31, 21]/255)
hold on
plot(tAxis(1:plotStep:end), pBubbleThermo(1:plotStep:end)/1e6, ...
     'k', 'LineWidth', 1)
hold off
xlabel ('$t$ [s]', 'Interpreter', 'latex')
ylabel ('$p$ [MPa]', 'Interpreter', 'latex')
legend({'Bubble inlet pressure', 'Bubble average pressure'}, ...
        'Interpreter', 'latex')

set(gca, 'XScale', 'log', 'YMinorTick', 'on')
global_format();
grid on
xlim([1e-4,1])
ylim([0, 5])

nexttile(tL, 4, [2,1])
xi = [sS.shuttle_position];
plot(tAxis(1:plotStep:end), xi, 'k', ...
    'LineWidth', 1);
xlabel ('$t$ [s]', 'Interpreter', 'latex')
ylabel ('$\xi$ [m]', 'Interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'YMinorTick', 'on')
global_format();
grid on
xlim([1e-4,1])
ylim([0, 0.06])

nexttile(tL, 6, [2,1])
% Plot bubble volume
% Use the decoupled continuation to plot bubble volume
tAxisContinuation = ...
    futuresBubbleModel(1).OutputArguments{1}.bubbleContinuationTime;
bubbleVol = 4/3*pi* ...
    futuresBubbleModel(1).OutputArguments{1}.bubbleContinuationState(1,:).^3;

plot(tAxisContinuation, bubbleVol, 'k', ...
    'LineWidth', 1);
xlabel ('$t$ [s]', 'Interpreter', 'latex')
ylabel ('$V$ [m${}^3$]', 'Interpreter', 'latex')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
global_format();
grid on
xlim([1e-4,1])
ylim([1e-5, 1e2])
set(gca, 'YTick', [1e-5, 1e0], 'YMinorTick', 'on')

nexttile(tL, 8, [2,1])
plot(tSample-tFirst, signals15/1e3, 'k', 'LineWidth', 1);
xlabel ('$t - r_1/c_\infty$ [s]', 'Interpreter', 'latex')
ylabel ('$\Delta p$ [kPa]', 'Interpreter', 'latex')

set(gca, 'XScale', 'log')
global_format();
grid on
xlim([1e-4,1])
ylim([-30, 70])
set(gca, 'YTick', [-30, 20, 70], 'YMinorTick', 'on')

set(gcf, 'Position', [674   338   595   609]);

%% Wrappers
function [p, tFirst] = sampleSignature(soln, metadata, tSample)
    % Returns p(t), and tFirst = r1/c_inf
    [pressureSignal, ~, ~, tFirst] = airgunSignature(soln, metadata);
    p = pressureSignal(tSample);
end

%% Compute at hydrophone
% Derived from airgunShuttleSignature file
function [pressureSignalFnTotal, tInterp, pressureSignal, tFirst] = ...
    airgunSignature(...
    solution, metadata, hydrophoneDepth, lateralSeparation)
if nargin < 3
    % Backward compatible default hydrophone depth
    hydrophoneDepth = 9;
    lateralSeparation = 6;
elseif nargin == 3
    error("Provide both hydrophone depth and lateral separation " + ...
          "(or neither for default).")
end

c_inf = 1482;     % Speed of sound in water [m/s]
rho_inf = 1000;   % Density of water [kg/m^3]
depth = metadata.paramAirgun.airgunDepth;

%% Define functions
RFn = @(t) bubbleRadiusFn(solution, t, size(solution.q,1));
RDotFn = @(t) bubbleVelocityFn(solution, t, size(solution.q,1));
RDotDotFn = @(t) bubbleAccelFn(solution, t, size(solution.q,1));
VDotFn = @(t) 4*pi*RFn(t).^2 .* RDotFn(t);
VDotDotFn = @(t) 4*pi*(...
    2*RFn(t) .* RDotFn(t).^2 ...
    + RFn(t).^2 .* RDotDotFn(t));

% Floor depth (ish)
floorDepth = 27;

% Using r as lateral distance
r1 = norm([lateralSeparation, depth-hydrophoneDepth]);
r2 = norm([lateralSeparation, depth+hydrophoneDepth]);% Ghost trig APPROX
r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepth]);% Hypothetical floor reflect

% Boolean value for using nonlinear, near-field pressure term
% from Keller and Kolodner
useNonlinearTerm = 0;

pressureDirect = @(t) rho_inf / (4*pi) * ...
    (VDotDotFn(t - r1/c_inf) / r1 ) - ...
    useNonlinearTerm * rho_inf / (32*pi^2*r1^4) * VDotFn(t - r1/c_inf).^2;
pressureGhost = @(t) rho_inf / (4*pi) * ...
    (VDotDotFn(t - r2/c_inf) / r2) - ...
    useNonlinearTerm * rho_inf / (32*pi^2*r2^4) * VDotFn(t - r2/c_inf).^2;
pressureThird = @(t) rho_inf / (4*pi) * ...
    (VDotDotFn(t - r3/c_inf) / r3) - ...
    useNonlinearTerm * rho_inf / (32*pi^2*r3^4) * VDotFn(t - r3/c_inf).^2;
pressureSignalFnTotal = @(t) pressureDirect(t) - pressureGhost(t) ...
    + pressureThird(t);

%% Compute pressure signal over a default grid
tInterp = linspace(metadata.tspan(1), metadata.tspan(end), 1000);
pressureSignal = pressureSignalFnTotal(tInterp);

%% Plot pressure signal if no output required
if nargout == 0
    plot(tInterp, pressureSignalFnTotal(tInterp));
end

tFirst = r1/c_inf;

%% [NEW] paper figures

[fullState, caseKeyContext] = ...
        airgunShuttlePostprocess( ...
        futuresBubbleModel(1).OutputArguments{1}, ...
        futuresBubbleModel(1).OutputArguments{2});

[agState, exception] = ...
    fullState(obj, q, t, bubble, shuttle, REVERT_MODEL, INCLUDE_ALL_PRIMITIVES)

end

function R = bubbleRadiusFn(solution, t, qSize)
    R = nan(size(t));
    
    for i = 1:length(t)
        if t(i) < 0
            % Initial quiescent signal
            R(i) = solution.bubbleContinuationState(1);
        elseif t(i) <= solution.soln.x(end)
            % Extract bubble info from coupled phase
            [sol, solDY] = deval(solution.soln, t(i));
            bubbleRIndex = qSize + 1;
            R(i) = sol(bubbleRIndex,:);
        else
            % Extract bubble info from uncoupled phase
            [sol, solDY] = deval(solution.solnBubbleContinuation, t(i));
            R(i) = sol(1,:);
        end
    end
end

function RDot = bubbleVelocityFn(solution, t, qSize)
    RDot = nan(size(t));
    
    for i = 1:length(t)
        if t(i) < 0
            % Initial quiescent signal
            RDot(i) = 0;
        elseif t(i) <= solution.soln.x(end)
            % Extract bubble info from coupled phase
            [sol, solDY] = deval(solution.soln, t(i));
            bubbleRDotIndex = qSize + 2;
            RDot(i) = sol(bubbleRDotIndex,:);
        else
            % Extract bubble info from uncoupled phase
            [sol, solDY] = deval(solution.solnBubbleContinuation, t(i));
            RDot(i) = sol(2,:);
        end
    end
end

function RDotDot = bubbleAccelFn(solution, t, qSize)
    % Initial quiescent signal
    RDotDot = nan(size(t));
    
    for i = 1:length(t)
        if t(i) < 0
            % Initial quiescent signal
            RDotDot(i) = 0;
        elseif t(i) <= solution.soln.x(end)
            % Extract bubble info from coupled phase
            [sol, solDY] = deval(solution.soln, t(i));
            bubbleRDotIndex = qSize + 2;
            RDotDot(i) = solDY(bubbleRDotIndex,:);
        else
            % Extract bubble info from uncoupled phase
            [sol, solDY] = deval(solution.solnBubbleContinuation, t(i));
            RDotDot(i) = solDY(2,:);
        end
    end
end