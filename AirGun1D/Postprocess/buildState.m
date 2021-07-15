% Builds state history struct

function exports = buildState( ...
    solution, metadata, fielddata)

%% Process input args
if nargin < 3
    dataAvailable = false;
else
    dataAvailable = true;
end

[fullState] = ...
        airgunShuttlePostprocess( ...
        solution, ...
        metadata);

exports = struct();

%% Compute pressure signature at receiver
tSample = linspace(0, 5, 15000);
signalModel = sampleSignature( ...
    solution, ...
    metadata, ...
    tSample);
exports.tSample = tSample;
exports.signalModel = signalModel;

%% Load data
% Set default values for field data x and y
timeDAQ = [];
signalData = [];

if dataAvailable
    dataDAQ = fielddata(25).entriesDAQ(4,1375:45000);
    timeDAQ = fielddata(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
    % Subtract dt to set first entry to t = 0
    timeDAQ = timeDAQ - timeDAQ(1);
    
    % Add static DAQ parameters
    DAQGain = 8;
    DAQSens = 1e5/7.6; % [Pa/V]
    signalData = DAQGain*DAQSens*dataDAQ;
    % Vertical shift signal data
    signalData = signalData - signalData(1);
    
    % Frequency domain information
    dwData = 1/(timeDAQ(end)-timeDAQ(1));
    omegaVecData = [0:1:(length(timeDAQ)/2-1), ...
        -length(timeDAQ)/2:1:-1] * dwData;
end

exports.timeDAQ = timeDAQ;
exports.signalData = signalData;

%% Extract mass flow rate
pS = [fullState.portStates];
exports.portMassFlow = [pS.massFlowPort];
exports.t = [fullState.t];

%% Frequency domain analysis

% Restrict tSample to timeDAQ if data is available
idxSampling = length(tSample);
if dataAvailable
    idxSampling = find(tSample >= timeDAQ(end), 1, 'first');
    idxSampling = 2*ceil(idxSampling/2);
end
tSampleFFT = tSample(1:idxSampling);
% Set up frequency axis
dw = 1/tSampleFFT(end);
omegaVec = [0:1:(length(tSampleFFT)/2-1), ...
    -length(tSampleFFT)/2:1:-1] * dw;

% Compute fft of model
exports.modelfft = fft(signalModel(1:idxSampling)) / length(omegaVec);
exports.omegaVecPositive = omegaVec(1:length(omegaVec)/2);
exports.modelfftPositive = exports.modelfft(1:length(omegaVec)/2);

% Set default values for fft of data
exports.datafft = [];
exports.omegaVecDataPositive = [];
exports.datafftPositive = [];
% Compute fft of data
if dataAvailable
    exports.datafft = fft(signalData) / length(omegaVecData);
    exports.omegaVecDataPositive = omegaVecData(1:length(omegaVecData)/2);
    exports.datafftPositive = exports.datafft(1:length(omegaVecData)/2);
end

%% End-pressure
eDS = [fullState.eulerDomainStates];
pAll = [eDS.p];
exports.p_LGrid = pAll(1,:);

% Prep experimental data (iNet)
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % Voltage baseline at 0 psi
% Manual data index input
begin_index = 99600; % (#25)
begin_index = 3700*5;
end_index = begin_index+1600*5;

% Extract iNet pressure data at wall boundary
% Time axis trimmed to first wall event
exports.iNetTimeAxis_p_L = [];
exports.iNetp_L = [];
if dataAvailable
    exports.iNetTimeAxis = fielddata(24).iNetTimeAxisP(begin_index:end_index) - ...
        fielddata(24).iNetTimeAxisP(begin_index);
    exports.iNetp_L = nominalminimum + ...
        psiPa_conversion*fielddata(24).iNetCh16Data(begin_index:end_index)';
end

%% Port outer pressure
exports.iNetTimeAxis_p_B = [];
exports.iNetp_B = [];
begin_index = 3660*5;
end_index = begin_index+1600*5;
if dataAvailable
    % Time axis trimmed to bubble event (earlier than wall event)
    exports.iNetTimeAxis_p_B = ...
        fielddata(24).iNetTimeAxisP(begin_index:end_index) ...
        - fielddata(24).iNetTimeAxisP(begin_index);
    exports.iNetp_B = ( ...
        psiPa_conversion*fielddata(24).iNetCh10Data(begin_index:end_index));
end

%% Bubble pressure
bS = [fullState.bubbleStates];
exports.bubblePressure = [bS.p];

end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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