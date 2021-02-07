%% Bubble and pressure signatures
% Wave propagation model

function [pressureSignalFn, tInterp, pressureSignal] = ...
    airgunShuttleSignature(solution, metadata)

r = 6;           % Distance from source to receiver [m]
c_inf = 1482;     % Speed of sound in water [m/s]
rho_inf = 1000;   % Density of water [kg/m^3]
depth = metadata.paramAirgun.airgunDepth;

%% Define functions
RFn = @(t) bubbleRadiusFn(solution.soln, t, size(solution.q,1));
RDotFn = @(t) bubbleVelocityFn(solution.soln, t, size(solution.q,1));
RDotDotFn = @(t) bubbleAccelFn(solution.soln, t, size(solution.q,1));

pressureDirect = @(t) rho_inf / (4*pi) * ...
    [RDotDotFn(t - r/c_inf) /r ];
pressureGhost = @(t) rho_inf / (4*pi) * ...
    [RDotDotFn(t - (r+2*depth)/c_inf) / (r+2*depth)];
pressureSignalFn = @(t) pressureDirect(t) - pressureGhost(t);

%% Compute pressure signal over a default grid
tInterp = linspace(metadata.tspan(1), metadata.tspan(end), 1000);
pressureSignal = pressureSignalFn(tInterp);

%% Plot pressure signal if no output required
if nargout == 0
    plot(tInterp, pressureSignalFn(tInterp));
end

end

%% Pressure data at closed end
function R = bubbleRadiusFn(odeSoln, t, qSize)
    % Initial quiescent signal
    R = nan(size(t));
    R(t < 0) = 0;

    [sol, solDY] = deval(odeSoln, t(t >= 0));
    bubbleRIndex = qSize + 1;
    R(t >= 0) = sol(bubbleRIndex,:);
end

function RDot = bubbleVelocityFn(odeSoln, t, qSize)
    % Initial quiescent signal
    RDot = nan(size(t));
    RDot(t < 0) = 0;
    
    [sol, solDY] = deval(odeSoln, t(t >= 0));
    bubbleRDotIndex = qSize + 2;
    RDot(t >= 0) = sol(bubbleRDotIndex,:);
end

function RDotDot = bubbleAccelFn(odeSoln, t, qSize)
    % Initial quiescent signal
    RDotDot = nan(size(t));
    RDotDot(t < 0) = 0;

    [sol, solDY] = deval(odeSoln, t(t >= 0));
    
    bubbleRDotIndex = qSize + 2;
    RDotDot(t >= 0) = solDY(bubbleRDotIndex,:);
end