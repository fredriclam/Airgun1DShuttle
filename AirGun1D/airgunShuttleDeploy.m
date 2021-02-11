% 
function [solution, metadata] = airgunShuttleDeploy(nx, coupleToShuttle)
if nargin == 1
    % Default: include shuttle model
    coupleToShuttle = true;
end

addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Simulation controls
% Simulation window [s]
% Suggested values:
% 0.100 to 0.600 s
tspan = [0; 0.200];
% Set flag for running shuttle-free model
runShuttleFreeFlag = false;

%% Set ambient parameters
airgunDepth = 10;                                 % Depth of airgun [m]
bubbleInitialVolume = 600;                        % Initial volume [cui]
airgunPressure = 1000;                            % Initial pressure in airgun [psi]

%% Set firing chamber parameters
airgunVolume = 20600;                             % Volume of airgun [cui]
airgunInnerDiameter = 10.0;                     % Inner diameter of airgun [in]
airgunCrossSecArea = pi*airgunInnerDiameter^2/4;  % Firing chamber cross-sectional area [in^2]
airgunLength = airgunVolume / ...
    airgunCrossSecArea * 0.0254;                  % Firing chamber length [m]
% Function prescribing firing chamber profile for future implementation
% Not used in current version.
airgunFiringChamberProfile = @(x) error(...
    'Not implemented. Placeholder for firing chamber profile function.');

%% Set mid chamber operation mode
% midChamberMode = 'limit-vented';
midChamberMode = 'limit-closed';

%% Set port parameters
airgunPortAreaRatio = 0.5;                        % Portion of lateral area covered by port [-]
airgunOuterDiameter = 11.2;                       % Outer diameter of firing chamber [in]
airgunPortLength = 2.5;                           % Length of port [in]
airgunPortArea = airgunPortAreaRatio * ...
    pi * airgunOuterDiameter * airgunPortLength;  % Effective port area [in^2]
% Shuttle parameters
shuttleBdryPenaltyStrength = 1e11;                % Linear elastic penalty term for shuttle [N/m]

%% Set operating chamber specifications
airCushionLength = 0.542*0.0254;                  % Length of closed air cushioning effect [m]
accelerationLength = (3.009-0.542)*0.0254;        % Length over which shuttle accelerates freely [m]

% Compression factor of air cushion as function of shuttle position
airgunOperatingChamberProfile = @(xi) (xi - accelerationLength < 0) * 1 ...
    + (xi - accelerationLength > 0) * ...
    (airCushionLength / (airCushionLength - (xi - accelerationLength)));

%% Housekeeping: Data packing
paramAirgun = struct(...
    'airgunPressure', airgunPressure, ...
    'airgunLength', airgunLength, ...
    'airgunCrossSecArea', airgunCrossSecArea, ...
    'airgunPortArea', airgunPortArea, ....
    'airgunDepth', airgunDepth, ...
    'airgunFiringChamberProfile', airgunFiringChamberProfile, ...
    'airgunOperatingChamberProfile', airgunOperatingChamberProfile, ...
    'bubbleInitialVolume', bubbleInitialVolume, ...
    'shuttleBdryPenaltyStrength', shuttleBdryPenaltyStrength, ...
    'midChamberMode', midChamberMode);
% Initialize metadata struct for documenting results
metadata = struct(...
    'paramAirgun', paramAirgun, ...
    'nx', nx, ...
    'tspan', tspan, ...
    'usingShuttleModel', coupleToShuttle);

%% Run solve for both models    
[solution, metadata] = ...
    runEulerCodeShuttleDual(nx, tspan, ...
                            paramAirgun, coupleToShuttle, metadata);