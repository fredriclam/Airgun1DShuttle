figure(16); clf;
tL = tiledlayout(5,3);

%% Header
% Add required dependencies from AirGun1D root.
% The following assumes the user is in the ./AirGun1D folder (otherwise,
% replace the following lines with the relative path to the respective
% folders.
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
futuresVar_mvar(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'bubbleModel', struct( ...
          'type', 'single', ...
          'M', 10, ...
          'alpha', 0.8), ...
        'extraOptions', struct( ...
          'shuttleAssemblyMass', 63/2 * .454) ...
));
futuresVar_mvar(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'bubbleModel', struct( ...
          'type', 'single', ...
          'M', 10, ...
          'alpha', 0.8), ...
        'extraOptions', struct( ...
          'shuttleAssemblyMass', 63*2 * .454) ...
));
futuresVar_mvar(3) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'bubbleModel', struct( ...
          'type', 'single', ...
          'M', 10, ...
          'alpha', 0.8), ...
        'extraOptions', struct( ...
          'shuttleAssemblyMass', 63/4 * .454) ...
));

masses = [7.15, 14.3, 28.6, 57.20];

%% Compute postprocess data (signals, fft of signal, wall pressure)
wait(futuresVar_mvar);
disp('Run complete.')
%% Memory management
solution_shlight = futuresVar_mvar(1).OutputArguments{1};
metadata_shlight = futuresVar_mvar(1).OutputArguments{2};
solution_shheavy = futuresVar_mvar(2).OutputArguments{1};
metadata_shheavy = futuresVar_mvar(2).OutputArguments{2};
solution_shfeath = futuresVar_mvar(3).OutputArguments{1};
metadata_shfeath = futuresVar_mvar(3).OutputArguments{2};
clear futuresVar_mvar;

%% Get full state history
[fullState_shlight, ~] = ...
        airgunShuttlePostprocess( ...
        solution_shlight, ...
        metadata_shlight);
[fullState_shheavy, ~] = ...
        airgunShuttlePostprocess( ...
        solution_shheavy, ...
        metadata_shheavy);
[fullState_shfeath, ~] = ...
        airgunShuttlePostprocess( ...
        solution_shfeath, ...
        metadata_shfeath);

%% Mass rate out
nexttile(tL, 7, [3,3]);
agtools.plotFiringChamber_exit(fullState_shfeath);
hold on
agtools.plotFiringChamber_exit(fullState_shlight);
agtools.plotFiringChamber_exit(fullState);
agtools.plotFiringChamber_exit(fullState_shheavy);
hold off

ch = get(gca,'Children');
set(ch(4), 'Color', 0.85*[1 1 1], 'LineStyle', '--', 'LineWidth', 1.5);
set(ch(3), 'Color', 0.6*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(2), 'Color', 0.3*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(1), 'Color', 0.0*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);

xlim([0, 25])
legend(arrayfun(@(m) sprintf("%.1f kg", m), masses), ...
       'Interpreter', 'tex', ...
       'location', 'best')
xlabel("{\it{t}} (ms)", "Interpreter", "tex")
ylim([0, 600])

set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
%% Shuttle position
temp.sS = [fullState_shfeath.shuttleStates];
shuttle_position_shfeath = [temp.sS.shuttle_position];
temp.sS = [fullState_shlight.shuttleStates];
shuttle_position_shlight = [temp.sS.shuttle_position];
temp.sS = [fullState_shheavy.shuttleStates];
shuttle_position_shheavy = [temp.sS.shuttle_position];
temp.sS = [fullState.shuttleStates];
shuttle_position_ref = [temp.sS.shuttle_position];

nexttile(tL, 2, [2,2]);
plot(1e3*[fullState_shfeath.t], shuttle_position_shfeath, 'r', ...
    'LineWidth', 1.5)
hold on
plot(1e3*[fullState_shlight.t], shuttle_position_shlight, 'g', ...
    'LineWidth', 1.5)
plot(1e3*[fullState.t], shuttle_position_ref, 'b', ...
    'LineWidth', 1.5)
plot(1e3*[fullState_shheavy.t], shuttle_position_shheavy, 'k', ...
    'LineWidth', 1.5)
hold off
xlim([0, 300])
ylim([0.00, 0.080])

ch = get(gca,'Children');

xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 14)
ylabel('$\xi$ (m)', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', "YTick", [0, 0.04, 0.08]);
   
ch = get(gca,'Children');
set(ch(4), 'Color', 0.85*[1 1 1], 'LineStyle', '--', 'LineWidth', 1.5);
set(ch(3), 'Color', 0.6*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(2), 'Color', 0.3*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(1), 'Color', 0.0*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);

nexttile(tL, 1, [2,1]);
box on
plot(1e3*[fullState_shfeath.t], shuttle_position_shfeath, 'r', ...
    'LineWidth', 1.5)
hold on
plot(1e3*[fullState_shlight.t], shuttle_position_shlight, 'g', ...
    'LineWidth', 1.5)
plot(1e3*[fullState.t], shuttle_position_ref, 'b', ...
    'LineWidth', 1.5)
plot(1e3*[fullState_shheavy.t], shuttle_position_shheavy, 'k', ...
    'LineWidth', 1.5)
hold off
hold off
xlim([0, 5])
ylim([0.00, 0.040])
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
ch = get(gca,'Children');
set(ch(4), 'Color', 0.85*[1 1 1], 'LineStyle', '--', 'LineWidth', 1.5);
set(ch(3), 'Color', 0.6*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(2), 'Color', 0.3*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
set(ch(1), 'Color', 0.0*[1 1 1], 'LineStyle', '-', 'LineWidth', 1.5);
xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 14)
ylabel('$\xi$ (m)', 'Interpreter', 'latex', 'FontSize', 14)