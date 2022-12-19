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

%% Global parameters
nx = 40;

%% Bubble comparison
% Background process or parallel (requires O(10 GB) per thread)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(1);
end

futuresFiner = parfeval(pool, @airgunShuttleDeploy, 2, ...
    60, true);
futuresCoarser = parfeval(pool, @airgunShuttleDeploy, 2, ...
    30, true);

%% Memory management
disp('Awaiting results')
wait(futuresFiner)
wait(futuresCoarser)
solution_finer = futuresFiner.OutputArguments{1};
metadata_finer = futuresFiner.OutputArguments{2};
clear futuresFiner;
%%
solution_coarser = futuresCoarser.OutputArguments{1};
metadata_coarser = futuresCoarser.OutputArguments{2};
clear futuresCoarser;

%% Load full state
[fullState_finer] = ...
        airgunShuttlePostprocess( ...
        solution_finer, ...
        metadata_finer);
    
[fullState_coarser] = ...
        airgunShuttlePostprocess( ...
        solution_coarser, ...
        metadata_coarser);

disp('Full state loaded')
return 

%% Plot wall states
figure(501);

colorCoarse = [252, 212, 98]/255
colorRef = 0*[1, 1, 1];
colorFine = [130, 213, 255]/255;

[wall_coarser.t, wall_coarser.p, ~] = agtools.plotFiringChamber_Wall( ...
    fullState_coarser, 'p', HiTestData);
hold on
[wall_ref.t, wall_ref.p, ~] = agtools.plotFiringChamber_Wall( ...
    fullState, 'p', HiTestData);
[wall_finer.t, wall_finer.p, ~] = agtools.plotFiringChamber_Wall( ...
    fullState_finer, 'p', HiTestData);
hold off

tiledlayout(2,1);
nexttile;
plot(wall_coarser.t*1e3, wall_coarser.p/1e6, '-', 'LineWidth', 2.5, ...
    'Color', colorCoarse)
hold on
plot(wall_ref.t*1e3, wall_ref.p/1e6, '-', 'LineWidth', 1.5, ...
    'Color', colorRef);
plot(wall_finer.t*1e3, wall_finer.p/1e6, ':', 'LineWidth', 2, ...
    'Color', colorFine)
hold off

legend(["$n_x = 30$", "$n_x = 40$", "$n_x = 60$"], ...
        'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
set(gcf, 'Position', [629   177   620   535]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel("$p$ [MPa]", ...
       'Interpreter', 'latex', 'FontSize', 14);

% Mass rate out plot
nexttile(2);
agtools.plotFiringChamber_exit(fullState_coarser);
hold on
agtools.plotFiringChamber_exit(fullState);
agtools.plotFiringChamber_exit(fullState_finer);
hold off

ch = get(gca,'Children');
set(ch(3), 'Color', colorCoarse, 'LineStyle', '-', 'LineWidth', 2.5);
set(ch(2), 'Color', colorRef,    'LineStyle', '-', 'LineWidth', 1.5);
set(ch(1), 'Color', colorFine,   'LineStyle', ':', 'LineWidth', 2.0);

legend(["$n_x = 30$", "$n_x = 40$", "$n_x = 60$"], ...
        'Interpreter', 'latex')

% Plot wall-pressure inset
inset_axes = axes('Position',[.38 .69 .30 .18]);
box on
plot(wall_coarser.t*1e3, wall_coarser.p/1e6, '-', 'LineWidth', 2.5, ...
    'Color', colorCoarse)
hold on
plot(wall_ref.t*1e3, wall_ref.p/1e6, '-', 'LineWidth', 1.5, ...
    'Color', colorRef);
plot(wall_finer.t*1e3, wall_finer.p/1e6, ':', 'LineWidth', 2, ...
    'Color', colorFine)
hold off
xlim([36, 46])
ylim([2.2, 2.8])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 2.2:0.2:2.8);

% Plot mass-out inset
inset_axes = axes('Position',[.38 .26 .30 .16]);
box on
agtools.plotFiringChamber_exit(fullState_coarser);
hold on
agtools.plotFiringChamber_exit(fullState);
agtools.plotFiringChamber_exit(fullState_finer);
xlabel ''
ylabel ''
hold off
xlim([55, 65])
ylim([120, 180])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');

ch = get(gca,'Children');
set(ch(1), 'Color', colorFine,   'LineStyle', ':', 'LineWidth', 2.0);
set(ch(2), 'Color', colorRef,    'LineStyle', '-', 'LineWidth', 1.5);
set(ch(3), 'Color', colorCoarse, 'LineStyle', '-', 'LineWidth', 2.5);