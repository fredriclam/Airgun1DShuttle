%% Runs the reference case in main Figures
% Generates firing chamber figures and so on

%% Header
assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
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

futuresReference = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', 'single', 'extraOptions', ...
    struct('TInitial', 288)));
futuresNoShuttle = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, false, ...
    struct('bubbleModel', 'single'));

%% Memory management
solution_reference = futuresReference.OutputArguments{1};
metadata_reference = futuresReference.OutputArguments{2};
clear futuresReference;
%%
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

%% Plot wall states
figure(1); clf;
subplot(3,1,1);
[t, p, pData] = agtools.plotFiringChamber_Wall(fullState, 'p', HiTestData);
subplot(3,1,2);
[t, T, TData] = agtools.plotFiringChamber_Wall(fullState, 'T', HiTestData);

subplot(3,1,3);
get_exp_entropy = @(p, T) ...
    T.^metadata_reference.discretization.physConst.gamma ./ ...
    p.^(metadata_reference.discretization.physConst.gamma-1);
es0 = get_exp_entropy(p(1,1), T(1,1));
get_entropy = @(p, T) log(get_exp_entropy(p,T) ./ es0);

plot(1e3*t, get_entropy(p, T) / es0, 'k--', 'LineWidth', 1);
hold on
plot(1e3*t, get_entropy(1e6*pchip(pData.t, pData.y, 1e3*t), ...
                        pchip(TData.t, TData.y, 1e3*t)), ...
     'b-', 'LineWidth', 1);
hold off

% Manual format
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel("$\gamma \ln \left( T/T_0 \right)" + ...
       "- (\gamma - 1) \ln \left( p/p_0 \right)$", ...
       'Interpreter', 'latex', 'FontSize', 11);
legend({'Model', 'Data'}, ...
      'Interpreter', 'latex', 'FontSize', 13);

%% Contour plot
figure(2); clf;
agtools.plotFiringChamber_xt( ...
fullState, metadata_reference, caseKeyContext);

%% Mass rate out plot
figure(3); clf;
agtools.plotFiringChamber_exit(fullState);

%% Bubble
figure(4); clf;
agtools.plotBubble(...
    fullState, ...
    solution_reference, ...
    metadata_reference);

%% Hydrophone signal and FFT
figure(5); clf;

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
    plot(timeDAQ(1:3000), ...
         dataTrack(i,1:3000), ...
         'LineWidth', 1);
    hold on
end
xlim([0, 15])
drawnow

disp("Data (acoustic pressure) is shown on figure.")

hold on
tSample = linspace(0,0.1,2000);
modelHydrophoneSignal = nan(4,length(tSample));
signalAcousticFns = {};
for i = 1:length(hydrophoneDepths)
    daqIdx = sortidx(i);
    [~, p, dists, funcs] = agtools.plotSignal(...
        fullState, ...
        solution_reference, ...
        metadata_reference, ...
        hydrophoneDepths(daqIdx), ...
        lateralSeparation, ...
        tSample);
    modelHydrophoneSignal(i,:) = p;
    signalAcousticFns{i} = @(t) dists.r1 * funcs.pFn(t)/1e5;
    disp("Plot model signal, step " + i)
end

xlim([0, 15])
legend(...
    [ ...
        arrayfun(@(x) ""+x+" m depth, data", hydrophoneDepths(sortidx)), ...
        arrayfun(@(x) ""+x+" m depth, model", hydrophoneDepths(sortidx))
    ], ...
    'Interpreter', 'latex', 'Location', 'northwest');
ylabel('$r_1 \Delta p$ [bar $\cdot$ m]', 'Interpreter', 'latex')
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);
hold off

% Set custom colours for plot
colScheme =  [0         0         0
    0.3137    0.3137    0.3137
    0.5020    0.5020    0.5020
    0.8000    0.8000    0.8000
    0.3490    0.1216    0.1216
    0.6588    0.1686    0.1686
    1.0000    0.0667    0.0667
    1.0000    0.4784    0.4784];
ch = get(gca,'children');
for i = 1:8
    set(ch(i), 'color', colScheme(i,:))
end

% Hydrophone FFT
figure(6); clf;

% do fft(t, p * r1), spectrum for signals spherically corrected to 1m
[wAxis, modelfft] = ...
    agtools.plotFFT((timeDAQ-timeDAQ(1))*1e-3, dataTrack(3,:)*1e5, ...
    @(t) 1e5*signalAcousticFns{3}(t));

xlabel("$f$ [Hz]", 'Interpreter', 'latex');
ylabel("SPL [dB re $1~\mu$Pa]", 'Interpreter', 'latex');
xlim([1e-1, 1e3]);
ylim([100, 220])
legend({'Model (9 m depth)', 'Data (9 m depth)'}, ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);
grid on

ch = get(gca,'children');
set(ch(1), 'color', [0.6588    0.1686    0.1686]);

%% Long window hydrophone plot
figure(91); clf;

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
    plot(timeDAQ(1:80000), ...
         dataTrack(i,1:80000), ...
         'LineWidth', 1);
    hold on
end
xlim([0, 2000])
drawnow

disp("Data (acoustic pressure) is shown on figure.")

hold on
tSample = linspace(0,2,500);
modelHydrophoneSignal = nan(4,length(tSample));
signalAcousticFns = {};
for i = 1:length(hydrophoneDepths)
    daqIdx = sortidx(i);
    [~, p, dists, funcs] = agtools.plotSignal(...
        fullState, ...
        solution_reference, ...
        metadata_reference, ...
        hydrophoneDepths(daqIdx), ...
        lateralSeparation, ...
        tSample);
    modelHydrophoneSignal(i,:) = p;
    signalAcousticFns{i} = @(t) dists.r1 * funcs.pFn(t)/1e5;
    disp("Plot model signal, step " + i)
end

legend(...
    [ ...
        arrayfun(@(x) ""+x+" m depth, data", hydrophoneDepths(sortidx)), ...
        arrayfun(@(x) ""+x+" m depth, model", hydrophoneDepths(sortidx))
    ], ...
    'Interpreter', 'latex', 'Location', 'best');
ylabel('$r_1 \Delta p$ [bar $\cdot$ m]', 'Interpreter', 'latex')
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);
hold off

% Set custom colours for plot
colScheme =  [0         0         0
    0.3137    0.3137    0.3137
    0.5020    0.5020    0.5020
    0.8000    0.8000    0.8000
    0.3490    0.1216    0.1216
    0.6588    0.1686    0.1686
    1.0000    0.0667    0.0667
    1.0000    0.4784    0.4784];
ch = get(gca,'children');
for i = 1:8
    set(ch(i), 'color', colScheme(i,:))
end

%% Compare bubble states for shuttle model  with instant-open model
figure(7); clf

tL = agtools.plotBubble_light( ...
    fullState, ...
    solution_reference, ...
    metadata_reference);
agtools.plotBubble_light( ...
    fullState_noshuttle, ...
    solution_noshuttle, ...
    metadata_noshuttle, ...
    tL, [48, 186, 69]/255);

nexttile(tL, 1)
legend({'Coupled model', 'Instant-open'}, ...
        'Interpreter', 'latex', 'location', 'northwest')

%% Compare wall states for shuttle model with instant-open model
figure(8); clf
subplot(2,1,1);
agtools.plotFiringChamber_Wall(fullState, 'p');
hold on
agtools.plotFiringChamber_Wall(fullState_noshuttle, 'p', HiTestData);
drawnow;
% Recolor the new momdel line
line_h = get(gca,'Children');
line_h(2).Color = [48, 186, 69]/255;

legend({'Coupled model', 'Instant-open', 'Data'}, ...
'Interpreter', 'latex', 'location', 'northeast')
%% Double contour-plot
figure(9); clf
[C2, h2] = agtools.plotFiringChamber_xt_nofill( ...
    fullState_noshuttle, metadata_noshuttle, caseKeyContext, ...
    [48, 186, 69]/255);
clabel(C2,h2,[6.5, 4, 2.0, 1.5], 'Color', [92, 219, 42]/255, ...
    'FontSize', 12, 'Interpreter', 'latex')
hold on
[C1, h1] = agtools.plotFiringChamber_xt_nofill( ...
    fullState, metadata_reference, caseKeyContext, ...
    [0, 0, 0]);
clabel(C1,h1,[6.5, 4.5, 4, 2.5, 2.0, 1.5], 'labelspacing', 160, ...
    'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 50])

legend([h1, h2], {'Coupled model', 'Instant-open model'}, ...
       'location', 'southwest', ...
       'Interpreter', 'latex', ...
       'FontSize', 13)