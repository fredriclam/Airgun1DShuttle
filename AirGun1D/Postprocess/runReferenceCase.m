%% Runs the reference case in main Figures
% Generates firing chamber figures and so on

%% Header
% Add required dependencies from AirGun1D root.
% The following assumes the user is in the ./AirGun1D folder (otherwise,
% replace the following lines with the relative path to the respective
% folders.
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

% Attempt to load field data
try
    load 'HiTestData_v1'
catch
    warning('Field trial data not found.')
end

%% Global parameters
nx = 40;

%% Bubble comparison
% Background process or parallel (requires O(10 GB) per thread)
% If the parallel toolbox is not available, airgunShuttleDeploy can be
% called directly with the shown arguments.
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(1);
end

futuresReference = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single', 'extraOptions', ...
    struct('TInitial', 288)));
futuresNoShuttle = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, false, ...
    struct('bubbleModel', 'single'));

%% Memory management
wait(futuresReference);
solution_reference = futuresReference.OutputArguments{1};
metadata_reference = futuresReference.OutputArguments{2};
clear futuresReference;
%%
wait(futuresNoShuttle);
solution_noshuttle = futuresNoShuttle.OutputArguments{1};
metadata_noshuttle = futuresNoShuttle.OutputArguments{2};
clear futuresNoShuttle;

%% Load full state
[fullState, caseKeyContext] = ...
        airgunShuttlePostprocess( ...
        solution_reference, ...
        metadata_reference);
[fullState_noshuttle] = ...
        airgunShuttlePostprocess( ...
        solution_noshuttle, ...
        metadata_noshuttle);