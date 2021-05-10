assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Global parameters
nx = 40;

%% Bubble comparison
% Background process or parallel (requires O(10 GB) per thread)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(2);
end

futuresBubbleModel(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single'));
futuresBubbleModel(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'quad'));

%% Design parameters
% Design matrix
portAreaRatioVector = 110/180 * ...
    [0.25, 0.8, 1, 1.2, 2.0, 4.0, 8.0];
pressureVector = [400, 600, 800, 1000, 2000, 3000]; % psi

for i = 1:length(portAreaRatioVector)
    portAreaRatio = portAreaRatioVector(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futuresPortAreaRatio(i) = parfeval(pool, @airgunShuttleDeploy, 2, ...
        nx, true, ...
        struct('airgunPortAreaRatio', portAreaRatio));
end
for i = 1:length(pressureVector)
    airgunPressure = pressureVector(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futuresPressures(i) = parfeval(pool, @airgunShuttleDeploy, 2, ...
        nx, true, ...
        struct('airgunPressure', airgunPressure));
end

%% Process full state for reference cases

[fullStateSingle, caseKeyContextSingle] = ...
        airgunShuttlePostprocess( ...
        futuresBubbleModel(1).OutputArguments{1}, ...
        futuresBubbleModel(1).OutputArguments{2});
[fullStateQuad, caseKeyContextQuad] = ...
        airgunShuttlePostprocess( ...
        futuresBubbleModel(2).OutputArguments{1}, ...
        futuresBubbleModel(2).OutputArguments{2});
    
%% Acoustic pressure signals
tSample = linspace(0, 5, 15000);

signalsSingle = sampleSignature( ...
        futuresBubbleModel(1).OutputArguments{1}, ...
        futuresBubbleModel(1).OutputArguments{2}, ...
        tSample);
signalsQuad = sampleSignature( ...
        futuresBubbleModel(2).OutputArguments{1}, ...
        futuresBubbleModel(2).OutputArguments{2}, ...
        tSample);
for i = 1:length(portAreaRatioVector)
    signalsPortAreaRatio{i} = sampleSignature( ...
        futuresPortAreaRatio(i).OutputArguments{1}, ...
        futuresPortAreaRatio(i).OutputArguments{2}, ...
        tSample);
end
for i = 1:length(pressureVector)
    signalsPressures{i} = sampleSignature( ...
        futuresPressures(i).OutputArguments{1}, ...
        futuresPressures(i).OutputArguments{2}, ...
        tSample);
end

%% Peaks (data)
dataDAQ = HiTestData(25).entriesDAQ(4,1375:45000);
% dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));

figure(201);
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
signalData = DAQGain*DAQSens*dataDAQ;
plot(timeDAQ, signalData, '.m');

dwData = 1/(timeDAQ(end)-timeDAQ(1));
omegaVecData = [0:1:(length(timeDAQ)/2-1), ...
    -length(timeDAQ)/2:1:-1] * dwData;

%%
hold on
plot(tSample, signalsSingle);
plot(tSample, signalsQuad);
hold off


%% Frequency domain analysis
dw = 1/tSample(end);
omegaVec = [0:1:(length(tSample)/2-1), ...
    -length(tSample)/2:1:-1] * dw;

figure(110); clf
loglog(omegaVecData, abs(fft(signalData)))
hold on
loglog(omegaVec, abs(fft(signalsSingle)))
loglog(omegaVec, abs(fft(2*signalsQuad)))
hold off
legend(["Data", "Single bubble", "Quad"])

figure(111); clf
for i = 1:length(portAreaRatioVector)
    loglog(omegaVec, abs(fft(signalsPortAreaRatio{i})))
    hold on
end
legend("portAreaRatio = " + portAreaRatioVector(1:end) + "")

figure(112); clf
for i = 1:length(pressureVector)
    loglog(omegaVec, abs(fft(signalsPressures{i})))
    hold on
end
legend("p_{init} = " + pressureVector(1:end) + " psi")

%% End-pressure
figure(301); clf;

% Prep PDE data
new.eDS = [fullStateSingle.eulerDomainStates];
new.pAll = [new.eDS.p];
new.p_L = new.pAll(1,:);
% Prep experimental data
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % Voltage baseline at 0 psi
% Manual data index input
begin_index = 3700*5;
end_index = begin_index+1600*5;

% Plot pressure [MPa] vs. t [ms]
plot(1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisP(begin_index)), ...
    1e-6* (nominalminimum + ...
    psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index))...
    , 'b-', 'LineWidth', 1);
hold on
plot(1000*[fullStateSingle.t], 1e-6*new.p_L, 'k-', 'LineWidth', 1);

hold off

xlim([0, 500]);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Data', 'Model $p_\mathrm{L}$', 'Model $p_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')


%% Wrappers
function p = sampleSignature(soln, metadata, tSample)
    pressureSignal = airgunSignature(soln, metadata);
    p = pressureSignal(tSample);
end

%% Compute at hydrophone
% Derived from airgunShuttleSignature file
function [pressureSignalFnTotal, tInterp, pressureSignal] = ...
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