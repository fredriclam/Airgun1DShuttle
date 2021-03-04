%% Bubble and pressure signatures
% Wave propagation model

function [pressureSignalFnTotal, tInterp, pressureSignal] = ...
    airgunShuttleSignature(...
    solution, metadata, hydrophoneDepth, lateralSeparation)

if nargin < 3
    % Backward compatible default hydrophone depth
    hydrophoneDepth = 9;
    lateralSeparation = 6;
elseif nargin < 4
    hydrophoneDepth = 9;
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
floorDepth = 28;

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