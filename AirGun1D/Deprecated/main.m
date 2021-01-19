% Main script for orchestrating code runs.
% 
% Adapted from WIP test_launch_script_v1.m script.

% Launch script for testing models with diff from [LW2019]

clear;
clc;
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Set parameters
nx = 100;       % Number of grid points per 1 m of air gun length

r = 10;                                           % Distance from source to receiver [m]
c_inf = 1482;                                     % Speed of sound in water [m/s]
rho_inf = 1000;                                   % Density of water [kg/m^3]
airgunPressure = 1000;                            % [psi]
% airgunLength = 2000 / (pi*10.020^2/4) * 0.0254; % [m] For ~0.6 m *****
airgunLength = 20600 / (pi*10.020^2/4) * 0.0254;  % [m] For ~0.6 m *****

% covers 0.5 of cyl
airgunPortArea = 0.5 * 2.5 * pi * 11.2;           % [in^2] % OD; covers some % of cyl
% airgunPortArea = 10 * 2.5 * 0.50 * pi * 11;     % [in^2] % OD;
airgunCrossSecArea = pi*10.020^2/4;               % [in^2]
airgunDepth = 10;                                 % [m]
bubbleInitialVolume = 600;                        % [cui]
shuttleBdryPenaltyStrength = 1e11; % [N/m]

% Function prescribing firing chamber profile. Not used in current version.
airgunFiringChamberProfile = @(x) error(...
    'Not implemented. Placeholder for firing chamber profile function.');

accelerationLength = (3.009-0.542)*0.0254;        % [m]
airCushionLength = 0.542*0.0254;                  % [m]
% Compression factor as function of shuttle position
airgunOperatingChamberProfile = @(xi) (xi - accelerationLength < 0) * 1 ...
    + (xi - accelerationLength > 0) * ...
    (airCushionLength / (airCushionLength - (xi - accelerationLength)));

% Run solve for both models
[sol, q, bubble, shuttle, plug, ...
    q2, bubble2, shuttle2, plug2, monitorStates] = runEulerCodeShuttleDual(nx, ...
    airgunPressure, airgunLength, airgunCrossSecArea, ...
    airgunPortArea, airgunDepth, airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength);
t = sol.x; % time

%% Postprocess 1: Plot closed chamber evolution
monitor_postprocessing(monitorStates)