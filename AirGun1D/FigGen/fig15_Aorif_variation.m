figure(15); clf;

%% Header
assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Runs
% Set parallel pool
% If parallel toolbox is not installed, comment out the following and
% replace parfeval with airgunShuttleDeploy( nx, [true|false], [...] ).
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(2);
end
nx = 40;

%% Coupled, with power-law pressure in bubble
futuresVar_orif(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'extraOptions', struct( ...
          'OpRearOrificeArea', 0) ...
));
futuresVar_orif(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'extraOptions', struct( ...
          'OpRearOrificeArea', 1e-5) ...
));
futuresVar_orif(3) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'extraOptions', struct( ...
          'OpRearOrificeArea', 1e-4) ...
));

orif_areas = [0, 1.37e-6, 1e-5, 1e-4];
%% Compute postprocess data (signals, fft of signal, wall pressure)
wait(futuresVar_orif);
disp('Run complete.')
%% Memory management
solution_orif0 = futuresVar_orif(1).OutputArguments{1};
metadata_orif0 = futuresVar_orif(1).OutputArguments{2};
solution_10micro = futuresVar_orif(2).OutputArguments{1};
metadata_10micro = futuresVar_orif(2).OutputArguments{2};
solution_100micro = futuresVar_orif(3).OutputArguments{1};
metadata_100micro = futuresVar_orif(3).OutputArguments{2};
clear futuresVar_orif;

%% Get full state history
[fullState_orif0, ~] = ...
        airgunShuttlePostprocess( ...
        solution_orif0, ...
        metadata_orif0);
[fullState_10micro, ~] = ...
        airgunShuttlePostprocess( ...
        solution_10micro, ...
        metadata_10micro);
[fullState_100micro, ~] = ...
        airgunShuttlePostprocess( ...
        solution_100micro, ...
        metadata_100micro);

%% Analysis
tL = tiledlayout(2,1);

%% Panel a

nexttile(tL,1);
plotShuttlePos(fullState_orif0);
hold on
plotShuttlePos(fullState);
plotShuttlePos(fullState_10micro);
plotShuttlePos(fullState_100micro);
hold on

set(gca, "FontSize", 15, "LineWidth", 1)
xlabel("{\it{t}} (ms)", "FontSize", 15)
ylabel("\xi (m)", "FontSize", 15)

%% Panel b

nexttile(tL,2);
plotP_L(fullState_orif0);
hold on
plotP_L(fullState);
plotP_L(fullState_10micro);
plotP_L(fullState_100micro);
hold on

set(gca, "FontSize", 15, "LineWidth", 1)
xlabel("{\it{t}} (ms)", "FontSize", 15)
ylabel("{\it{p}}_{L} (MPa)", "FontSize", 15)

if exist("HiTestData", "var")
    hold on
    % Prep field data
    psiPa_conversion = 6894.75729;
    % Voltage baseline at 0 psi
    nominalminimum = 1000*psiPa_conversion; 
    % Manual data index input
    begin_index = 3700*5;
    end_index = begin_index+1600*5;
    fieldData = struct();
    fieldData.t = 1e3 * ...
        (HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
         HiTestData(24).iNetTimeAxisP(begin_index));
    fieldData.y = 1e-6* (nominalminimum + ...
        psiPa_conversion * ...
        HiTestData(24).iNetCh16Data(begin_index:end_index));
    plot(fieldData.t, fieldData.y);
    hold off
    
    legend(["0 mm^2", "1.37 mm^2", "10 mm^2", "100 mm^2", "Data"])
else
    legend(["0 mm^2", "1.37 mm^2", "10 mm^2", "100 mm^2"])
end
xlim([0, 300]);

function [t, xi] = plotShuttlePos(fullState)
sS = [fullState.shuttleStates];
xi = [sS.shuttle_position];
t = [fullState.t];
plot(1e3*t,xi,"LineWidth",1);
end

function [t, p_L] = plotP_L(fullState)
eDS = [fullState.eulerDomainStates];
p = [eDS.p];
p_L = p(1,:);
t = [fullState.t];
plot(1e3*t,p_L/1e6,"LineWidth",1);
end