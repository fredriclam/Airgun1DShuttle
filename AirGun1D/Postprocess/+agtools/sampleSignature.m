% Samples the hydrophone signal
%
% Outputs:
%   p: pressure signal interpolated at tSample
%   funcs: struct with:
%     VFn: Function handle for V = V(t)
%     VDotFn: Function handle for VDot = (d/dt) V(t)
%     VDotDotFn: Function handle for VDot = (d/dt) (d/dt) V(t)

function [p, funcs] = sampleSignature( ...
    soln, metadata, tSample, hydrophoneDepth, lateralSeparation)
    if nargin <= 2
        tSample = linspace(metadata.tspan(1), metadata.tspan(end), 1000);
    end
    
    if nargin <= 3
        % Backward compatible default hydrophone depth
        hydrophoneDepth = 9;
        lateralSeparation = 6;
    elseif nargin == 4
        error("Provide both hydrophone depth and lateral separation " + ...
              "(or neither for default).")
    end

    [pressureSignal, funcs] = airgunSignature( ...
        soln, metadata, hydrophoneDepth, lateralSeparation);
    p = pressureSignal(tSample);
    
    %% Plot pressure signal if no output required
    if nargout == 0
        plot(tSample, p);
    end
end

%% Compute at hydrophone
% Derived from airgunShuttleSignature file
function [pressureSignalFnTotal, funcs] = ...
    airgunSignature(...
    solution, metadata, hydrophoneDepth, lateralSeparation)

c_inf = 1482;     % Speed of sound in water [m/s]
rho_inf = 1000;   % Density of water [kg/m^3]
depth = metadata.paramAirgun.airgunDepth;

%% Define functions
RFn = @(t) bubbleRadiusFn(solution, t, size(solution.q,1));
RDotFn = @(t) bubbleVelocityFn(solution, t, size(solution.q,1));
RDotDotFn = @(t) bubbleAccelFn(solution, t, size(solution.q,1));
VFn = @(t) 4/3*pi*RFn(t).^3;
VDotFn = @(t) 4*pi*RFn(t).^2 .* RDotFn(t);
VDotDotFn = @(t) 4*pi*(...
    2*RFn(t) .* RDotFn(t).^2 ...
    + RFn(t).^2 .* RDotDotFn(t));

% Estimate of effective floor depth for reflection
floorDepth = 27;

% Using r as lateral distance
r1 = norm([lateralSeparation, depth-hydrophoneDepth]);
r2 = norm([lateralSeparation, depth+hydrophoneDepth]);%
r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepth]);

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

funcs = struct( ...
    'VFn', VFn, ...
    'VDotFn', VDotFn, ...
    'VDotDotFn', VDotDotFn, ...
    'pressureDirect', pressureDirect, ...
    'pFn', pressureSignalFnTotal ...
);
end

function R = bubbleRadiusFn(solution, t, qSize)
    R = nan(size(t));
    
    for i = 1:length(t)
        if t(i) < solution.soln.x(1)
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
        if t(i) < solution.soln.x(1)
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
        if t(i) < solution.soln.x(1)
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