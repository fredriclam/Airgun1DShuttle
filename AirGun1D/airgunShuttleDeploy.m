% Wrapper for running airgun-shuttle coupled model. 
% Sets default parameters, and overrides with valid fields provided through
% `options` by unpacking the variables into the local namespace.
% 
% 
function [solution, metadata] = airgunShuttleDeploy(nx, coupleToShuttle, options)
if nargin == 1
    % Default: include shuttle model
    coupleToShuttle = true;
end
if nargin >= 3
    useOverrideOptions = true;
else
    useOverrideOptions = false;
end

%% Simulation controls
% Simulation window [s]
% Suggested values:
% [0, and 0.100 to 0.600 s]
if nargin >= 3 && isfield(options, "leakTime")
    tspan = [-options.leakTime, 0.300];
else
    tspan = [0; 0.300];
end
% Set flag for running shuttle-free model
runShuttleFreeFlag = false;

%% Set bubble model
bubbleModel = 'single';

%% Set ambient parameters
airgunDepth = 10;                                 % Depth of airgun [m]
bubbleInitialVolume = 10;                         % Initial volume [cui] -- changed
if ~coupleToShuttle
    % Pick larger seed bubble volume to control initial bubble pressure
    % increase and maintain backward compatibility with Watson et al. 2019
    bubbleInitialVolume = 600;
end
airgunPressure = 1000;                            % Initial pressure in airgun [psi]

%% Set firing chamber parameters
airgunVolume = 20600;                             % Volume of airgun [cui]
airgunInnerDiameter = 10.0;                       % Inner diameter of airgun [in]
% Function prescribing firing chamber profile for future implementation
% Not used in current version.
airgunFiringChamberProfile = @(x) error(...
    'Not implemented. Placeholder for firing chamber profile function.');
% Set mid chamber operation mode
% midChamberMode = 'limit-vented';
midChamberMode = 'limit-closed';
% Set linear elastic penalty term for solid-solid collision shuttle [N/m]
shuttleBdryPenaltyStrength = 1e11; 
% Deprecated parameters (constrained by design measurement)
airgunPortLength = NaN;

%% Set legacy operating chamber specs
airgunOperatingChamberProfile = @() error('Legacy argument used');
%% Extra options struct for passing through to airgunConfig
extraOptions = struct();

%% Evaluate options
% Unpacks all the provided options into the local scope, replacing
% and above presets.
if useOverrideOptions
    optionsKeys = fields(options);
    for i = 1:length(optionsKeys)
        key = optionsKeys{i};
        eval(sprintf('%s = options.(key);', key))
    end
end

%% Compute dependent parameters
airgunCrossSecArea = pi*airgunInnerDiameter^2/4;            % Firing chamber cross-sectional area [in^2]
airgunLength = airgunVolume / airgunCrossSecArea * 0.0254;  % Firing chamber length [m]
% Set port area from measurement [in^2]
airgunPortArea = 90;

%% Housekeeping: Data packing
paramAirgun = struct(...
    'airgunPressure', airgunPressure, ...
    'airgunLength', airgunLength, ...
    'airgunCrossSecAreaSqInch', airgunCrossSecArea, ...
    'airgunPortAreaSqInch', airgunPortArea, ....
    'airgunDepth', airgunDepth, ...
    'airgunFiringChamberProfile', airgunFiringChamberProfile, ...
    'airgunOperatingChamberProfile', airgunOperatingChamberProfile, ...
    'bubbleInitialVolume', bubbleInitialVolume, ...
    'shuttleBdryPenaltyStrength', shuttleBdryPenaltyStrength, ...
    'midChamberMode', midChamberMode, ...
    'airgunPortLength', airgunPortLength, ...
    'bubbleModel', bubbleModel);
% Initialize metadata struct for documenting numerical solution
metadata = struct(...
    'paramAirgun', paramAirgun, ...
    'nx', nx, ...
    'tspan', tspan, ...
    'usingShuttleModel', coupleToShuttle, ...
    'extraOptions', extraOptions);

%% Solve PDE
[solution, metadata] = ...
    runEulerCodeShuttleDual(nx, tspan, ...
                            paramAirgun, coupleToShuttle, metadata);