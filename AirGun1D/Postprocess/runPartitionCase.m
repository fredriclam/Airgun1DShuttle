% Runs several test cases for the 
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
    pool = parpool(1);
end

nx = 40;
%% Coupled, with power-law pressure in bubble
futures_partition(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'partition')));
%% Compute postprocess data (signals, fft of signal, wall pressure)
wait(futures_partition);
disp('Run complete.')
%% Memory management
solution_partition = futures_partition.OutputArguments{1};
metadata_partition = futures_partition.OutputArguments{2};
clear futures_partition;

%% Get full state history
[fullState_partition, caseKeyContext_partition] = ...
        airgunShuttlePostprocess( ...
        solution_partition, ...
        metadata_partition);

%% Hydrophone signal and FFT
figure(401); clf;

% Plot hydrophone data explicitly, multiplying by r1 (incident wave travel
% distance) to scale data near the peak pressure
floorDepth = 27;
lateralSeparation = 6;
depth = metadata_reference.paramAirgun.airgunDepth;
hydrophoneDepths = [6,12,3,9];
[~, sortidx] = sort(hydrophoneDepths);
tSignal = HiTestData(25).headerDAQ.SamplingInterval * ...
    ((1:HiTestData(25).headerDAQ.SampleCount)-1);
signalsData = HiTestData(25).entriesDAQ;
start_index = 1470;
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
manual_tshift = 3.8;
timeDAQ = 1e3*(tSignal(start_index:end)-tSignal(start_index)) ...
          + manual_tshift;
clear dataTrack;

for i = 1:4
    daqIdx = sortidx(i);
    r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
    r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
    r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
    
    dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
    plot(timeDAQ(1:12000), ...
         dataTrack(i,1:12000), ...
         'LineWidth', 1);
    hold on
end
xlim([0, 100])
drawnow

disp("Data (acoustic pressure) is shown on figure.")

hold on
tSample = linspace(0,0.3,6000);
modelHydrophoneSignal_partition = nan(4,length(tSample));
signalAcousticFns = {};
for i = 1:length(hydrophoneDepths)
    daqIdx = sortidx(i);
    [~, p, dists, funcs] = agtools.plotSignal(...
        fullState_partition, ...
        solution_partition, ...
        metadata_partition, ...
        hydrophoneDepths(daqIdx), ...
        lateralSeparation, ...
        tSample, [0,0,0], '-');
    drawnow
    modelHydrophoneSignal_partition(i,:) = p;
    signalAcousticFns{i} = @(t) dists.r1 * funcs.pFn(t)/1e5;
    disp("Plot model signal, step " + i)
end

legend(...
    [ ...
        arrayfun(@(x) ""+x+" m depth, data", hydrophoneDepths(sortidx)), ...
        arrayfun(@(x) ""+x+" m depth, partition model", hydrophoneDepths(sortidx))
    ], ...
    'Interpreter', 'latex', 'Location', 'northeast');
ylabel('$r_1 \Delta p$ [bar $\cdot$ m]', 'Interpreter', 'latex')
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);
hold off

% Set custom colours for plot
colScheme =  [0.3686    0.0902    0.3647
    0.5608    0.1333    0.5529
    0.8588    0.2078    0.8510
    0.9608    0.5294    0.9529
    0.3490    0.1216    0.1216
    0.6588    0.1686    0.1686
    1.0000    0.0667    0.0667
    1.0000    0.4784    0.4784];
ch = get(gca,'children');
for i = 1:8
    set(ch(i), 'color', colScheme(i,:), 'LineWidth', 1.5)
end