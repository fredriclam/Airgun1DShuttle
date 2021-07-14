%% Performs matrix study (main thread version)

%% Header
assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Serial firing
nx = 40;
portAreaRatioVector = 110/180 * [0.25, 0.8, 1, 1.2, 2.0, 4.0];
pressureVector = [400, 600, 800, 1000, 2000, 3000]; % psi

for i = 1:length(portAreaRatioVector)
    for j = 1:length(pressureVector)
        portAreaRatio = portAreaRatioVector(i);
        airgunPressure = pressureVector(j);
        
        % Call function [solution, metadata] = airgunShuttleDeploy(nx, ...)
        [solutionTemp, metadataTemp] = airgunShuttleDeploy(nx, true, ...
            struct('airgunPortAreaRatio', portAreaRatio, ...
                   'airgunPressure', airgunPressure, ...
                   'bubbleModel', 'single'));
        % Strip data
        dataMatrix{i,j} = dataReduction(solutionTemp);
        metadataMatrix{i,j} = metadataTemp;
        % Clear mem
        clear solutionTemp;
        clear metadataTemp;
    end
end

%% Analysis

%% Precompute hydrophone pressure signal
hydrophoneDepth = 6;
lateralSeparation = 10;
RTE = nan(6,6);

for i = 1:6
    for j = 1:6
        [pressureSignalFnTotal, ~, ~] = ...
            airgunSignatureReduced(...
            dataMatrix{i,j}, metadataMatrix{i,j}, ...
            hydrophoneDepth, lateralSeparation);
        %% Sampling
        % Set 0.2 second window for sampling
        % Linear interpolation of Rdotdot is a poor way to recover information
        % Need the ode23 soln object and deval strategy
        tSampleFine = linspace(0, 0.2, 4000);
        % tSampleCoarse = linspace(0.2, 5, 5000);
        tSampleCoarse = [];
        tSample = [tSampleFine, tSampleCoarse(2:end)];
        p = pressureSignalFnTotal(tSample);
        
        %%
        figure(801); clf;
        plot(tSample*1e3, p);
        drawnow
        xlim([0, 200])
        %% Rise time estimation
        RTE(i,j) = ...
            dataMatrix{i,j}.bubbleContinuationTime(find(diff(p) < 0, 1, 'first')) ...
            - dataMatrix{i,j}.bubbleContinuationTime(find(diff(p) > 0, 1, 'first'))
    end
end

%% Generate colour plot
psi2MPa = 6894.76e-6;
[meshx, meshy] = meshgrid(pressureVector*psi2MPa, ...
    portAreaRatioVector);
contourf(meshx, meshy, RTE*1e3, 'LineStyle', 'none')

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$p_0$ [MPa]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$A_\mathrm{port,max}/A_\mathrm{cs}$';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$t_\mathrm{r}$ [ms]';

%% Compute signature at hydrophone, for reduced data structures only
function [pressureSignalFnTotal, tInterp, pressureSignal] = ...
    airgunSignatureReduced(...
    solution, metadata, hydrophoneDepth, lateralSeparation)
c_inf = 1482;     % Speed of sound in water [m/s]
rho_inf = 1000;   % Density of water [kg/m^3]
depth = metadata.paramAirgun.airgunDepth;

%% Define functions

function R = bubbleRadiusFn(solution, t)
    R = nan(size(t));
    for i = 1:length(t)
        if t(i) < 0
            % Initial quiescent signal
            R(i) = solution.bubbleContinuationState(1);
        else
            % Extract bubble info from continuation (interp linearly)
            R(i) = interp1(solution.bubbleContinuationTime, ...
                           solution.bubbleContinuationState(1,:), ...
                           t(i));
        end
    end
end

function RDot = bubbleVelocityFn(solution, t)
    RDot = nan(size(t));
    for i = 1:length(t)
        if t(i) < 0
            % Initial quiescent signal
            RDot(i) = 0;
        else
            % Extract bubble info from continuation (interp linearly)
            RDot(i) = interp1(solution.bubbleContinuationTime, ...
                              solution.bubbleContinuationState(2,:), ...
                              t(i));
        end
    end
end

function RDotDot = bubbleAccelFn(solution, t)
    % Initial quiescent signal
    RDotDot = nan(size(t));
    
    for i = 1:length(t)
        if t(i) <= 0
            % Initial quiescent signal
            RDotDot(i) = 0;
        else
            ind1 = find(solution.bubbleContinuationTime < t(i), ...
                1, 'last');
            ind2 = find(solution.bubbleContinuationTime > t(i), ...
                1, 'first');
            assert(ind2 == ind1+1);
            dRdot = solution.bubbleContinuationState(2,ind2) ...
                  - solution.bubbleContinuationState(2,ind1);
            dt = solution.bubbleContinuationTime(ind2) ...
               - solution.bubbleContinuationTime(ind1);
            RDotDot(i) = dRdot/dt;
        end
    end
end

RFn = @(t) bubbleRadiusFn(solution, t);
RDotFn = @(t) bubbleVelocityFn(solution, t);
RDotDotFn = @(t) bubbleAccelFn(solution, t);
VDotFn = @(t) 4*pi*RFn(t).^2 .* RDotFn(t);
VDotDotFn = @(t) 4*pi*(...
    2*RFn(t) .* RDotFn(t).^2 ...
    + RFn(t).^2 .* RDotDotFn(t));



% Basin floor depth (estimate)
floorDepth = 27;

% Using r as lateral distance
r1 = norm([lateralSeparation, depth-hydrophoneDepth]);
r2 = norm([lateralSeparation, depth+hydrophoneDepth]);
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

%% Compute pressure signal over a default grid
tInterp = linspace(metadata.tspan(1), metadata.tspan(end), 1000);
pressureSignal = pressureSignalFnTotal(tInterp);

%% Plot pressure signal if no output required
if nargout == 0
    plot(tInterp, pressureSignalFnTotal(tInterp));
end

end