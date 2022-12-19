% Run the leak case to test hypothesis of seed bubble influencing the
% acoustic peak pressure overprediction.

nx = 40;
p_ksi = 6.895e+6;
[solution_leak_rr0, metadata_leak_rr0] = airgunShuttleDeploy(nx, true, ...
    struct(...
        'bubbleInitialVolume', 10, ...
        'bubbleModel', ...
            struct('type', 'single', 'waterProperties', struct(), 'flowreductionfactor', 0.6), ...
    'extraOptions', struct('TInitial', 288, 'pb', 1*198000), ... % consider also 'R',0.3,
    'leakTime', 0.0));
%%
[solution_leak_rr1, metadata_leak_rr1] = airgunShuttleDeploy(nx, true, ...
    struct(...
        'bubbleModel', ...
            struct('type', 'single', 'waterProperties', struct(), 'flowreductionfactor', 0.6), ...
    'extraOptions', struct('TInitial', 288, 'pb', 1*198000, 'R', 0.1), ...
    'leakTime', 0.0));
[solution_leak_rr2, metadata_leak_rr2] = airgunShuttleDeploy(nx, true, ...
    struct(...
        'bubbleModel', ...
            struct('type', 'single', 'waterProperties', struct(), 'flowreductionfactor', 0.6), ...
    'extraOptions', struct('TInitial', 288, 'pb', 1*198000, 'R', 0.2), ...
    'leakTime', 0.0));
[solution_leak_rr3, metadata_leak_rr3] = airgunShuttleDeploy(nx, true, ...
    struct(...
        'bubbleModel', ...
            struct('type', 'single', 'waterProperties', struct(), 'flowreductionfactor', 0.6), ...
    'extraOptions', struct('TInitial', 288, 'pb', 1*198000, 'R', 0.3), ...
    'leakTime', 0.0));
%%
[solution_leak_rr4, metadata_leak_rr4] = airgunShuttleDeploy(nx, true, ...
    struct(...
        'bubbleModel', ...
            struct('type', 'partition', 'waterProperties', struct(), 'flowreductionfactor', 1.0), ...
    'extraOptions', struct('TInitial', 288, 'pb', 1*198000, 'R', 0.2), ...
    'leakTime', 0.0));

%% Postprocess the solution

[fullState_leak_0] = ...
        airgunShuttlePostprocess( ...
        solution_leak_rr0, ...
        metadata_leak_rr0);    
[fullState_leak_3] = ...
        airgunShuttlePostprocess( ...
        solution_leak_rr3, ...
        metadata_leak_rr3);    
%%
[fullState_leak_4] = ...
        airgunShuttlePostprocess( ...
        solution_leak_rr4, ...
        metadata_leak_rr4);

%% Plot port pressure of leak vs. no leak (assume refilling)

ps0 = [fullState_leak_0.portStates];
ps3 = [fullState_leak_3.portStates];
bs0 = [fullState_leak_0.bubbleStates];
bs3 = [fullState_leak_3.bubbleStates];

figure(4); clf
subplot(1,2,1);
plot([fullState_leak_0.t], [ps0.pPort])
hold on
plot([fullState_leak_3.t], [ps3.pPort])
xlim([0,0.3])

subplot(1,2,2);
plot([fullState_leak_0.t], [bs0.R])
hold on
plot([fullState_leak_3.t], [bs3.R])
xlim([0,0.3])

%% Postprocess for acoustic signature
[p_aco_leak_rr0, funcs_leak_rr0] = agtools.sampleSignature( ...
    solution_leak_rr0, metadata_leak_rr0);
%%
[p_aco_leak_rr1, funcs_leak_rr1] = agtools.sampleSignature( ...
    solution_leak_rr1, metadata_leak_rr1);
[p_aco_leak_rr2, funcs_leak_rr2] = agtools.sampleSignature( ...
    solution_leak_rr2, metadata_leak_rr2);
[p_aco_leak_rr3, funcs_leak_rr3] = agtools.sampleSignature( ...
    solution_leak_rr3, metadata_leak_rr3);
[p_aco_leak_rr4, funcs_leak_rr4] = agtools.sampleSignature( ...
    solution_leak_rr4, metadata_leak_rr4);

%% Plot using arbitrary time window
gen_tsp = @(md) linspace(md.tspan(1), md.tspan(end), 1000);
tsp_rr0 = gen_tsp(metadata_leak_rr0);
tsp_rr1 = gen_tsp(metadata_leak_rr1);
tsp_rr2 = gen_tsp(metadata_leak_rr2);
tsp_rr3 = gen_tsp(metadata_leak_rr3);
tsp_rr4 = gen_tsp(metadata_leak_rr4);

p_out_rr0 = funcs_leak_rr0.pFn(tsp_rr0);
p_out_rr1 = funcs_leak_rr1.pFn(tsp_rr1);
p_out_rr2 = funcs_leak_rr2.pFn(tsp_rr2);
p_out_rr3 = funcs_leak_rr3.pFn(tsp_rr3);
p_out_rr4 = funcs_leak_rr4.pFn(tsp_rr4);

figure(13); clf;
plot(tsp_rr0, p_out_rr0);
hold on
plot(tsp_rr1, p_out_rr1);
plot(tsp_rr2, p_out_rr2);
plot(tsp_rr3, p_out_rr3);
plot(tsp_rr4, p_out_rr4);
hold off

return

%% Get longer signal for tuned model
tsp_tuned = linspace(metadata_leak_rr4.tspan(1), 10, 5000);
[p_out_tuned, funcs_tuned] = agtools.sampleSignature( ...
    solution_leak_rr4, metadata_leak_rr4, tsp_tuned);
figure(23); clf;
plot(tsp_tuned, p_out_tuned);

%% Plot with data
figure(6); clf;
floorDepth = 27;
lateralSeparation = 6;
depth = metadata_leak_rr0.paramAirgun.airgunDepth;
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

r1_96 = norm([6, depth-9]);
hold on
plot(tsp_rr0, p_out_rr0/1e5*r1_96, 'r');
plot(tsp_rr1, p_out_rr1/1e5*r1_96, 'g');
plot(tsp_rr2, p_out_rr2/1e5*r1_96, 'b');
plot(tsp_rr3, p_out_rr3/1e5*r1_96, 'm');
plot(tsp_rr4, p_out_rr4/1e5*r1_96, 'c');

% Plot DAQ
for i = 1:4
    daqIdx = sortidx(i);
    r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
    r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
    r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
    
    dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
    plot(1e-3*timeDAQ(1:12000), ...
         dataTrack(i,1:12000), ...
         'LineWidth', 1);
    hold on
end

xlim([0, 1.0])
drawnow

% Get dynamic initial condition
R0 = pchip(solution_leak_rr0.soln.x, solution_leak_rr0.bubble(1,:),0);
Rprime0 = pchip(solution_leak_rr0.soln.x, solution_leak_rr0.bubble(2,:),0);

get_IC_label = @(soln) "R_0 = " + pchip(soln.soln.x, soln.bubble(1,:),0) ...
    + ", R'_0 = " + pchip(soln.soln.x, soln.bubble(2,:),0);
legend( ...
    get_IC_label(solution_leak_rr0), ...
    get_IC_label(solution_leak_rr1), ...
    get_IC_label(solution_leak_rr2), ...
    get_IC_label(solution_leak_rr3) ...
)

%% Plot hydrophones against tuned model
figure(71); clf;

% Plot DAQ
for i = 1:4
    daqIdx = sortidx(i);
    r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
    r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
    r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
    
    dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
    plot(1e3*1e-3*timeDAQ(1:72000), ...
         dataTrack(i,1:72000), ...
         'LineWidth', 1);
    hold on
end
plot(1e3*tsp_tuned, p_out_tuned/1e5*r1_96);
xlim([0, 400])
xlabel("t (ms)")
ylabel("r_1 \Delta p (bar \cdot m)")
set(gca, "FontSize", 14, "LineWidth", 1.0, "XMinorTick", "on", "YMinorTick", "on")

%% Plot internal thermodynamics
figure(72); clf;

% Plot DAQ
V_opt = [solution_leak_rr4.bubbleContinuationState(1,:)].^3 * (4*pi/3);
E_opt = [solution_leak_rr4.bubbleContinuationState(4,:)];
pb_opt = (1.4 - 1) * E_opt ./ V_opt;
t_opt = [solution_leak_rr4.bubbleContinuationTime];
plot(1e3*t_opt, 1e-6*pb_opt);
xlim([0,15])
hold on

%% Overlay internal thermo data
% Prep field data
psiPa_conversion = 6894.75729;
% Voltage baseline at 0 psi
nominalminimum = 1e5 + 9.8*10*1000;
begin_index = 3700*5;
end_index = begin_index+1600*5;
fieldData = struct();
fieldData.t = 1e3 * ...
    (HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
     HiTestData(24).iNetTimeAxisP(begin_index));
fieldData.y = 1e-6* (nominalminimum + ...
    psiPa_conversion * ...
    HiTestData(24).iNetCh10Data(begin_index:end_index));
hold on
plot(fieldData.t, fieldData.y, ...
     'LineWidth', 1.5, ...
     'Color', [0.6588    0.1686    0.1686])
hold off
xlim([0, 15])

ylabel('{\it{p}} (MPa)', 'Interpreter', 'tex', 'FontSize', 15)
xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 15)
set(gca, 'FontSize', 15, ...
    'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);


%% Extended bubble radius
t_leak = [solution_leak.bubbleContinuationTime];
Rb_leak = [solution_leak.bubbleContinuationState(1,:)];
Vb_leak = 4*pi/3* [solution_leak.bubbleContinuationState(1,:)].^3;
t_ref = [solution_reference.bubbleContinuationTime];
Rb_ref= [solution_reference.bubbleContinuationState(1,:)];
Vb_ref = 4*pi/3* [solution_reference.bubbleContinuationState(1,:)].^3;


%% Plot using arbitrary time window
tsp_ref = linspace(metadata_reference.tspan(1), metadata_reference.tspan(end), 1000);
tsp1 = linspace(metadata_leak.tspan(1), metadata_leak.tspan(end), 1000);
tsp2 = linspace(metadata_leak2.tspan(1), metadata_leak2.tspan(end), 1000);

p_outref = funcs_ref.pFn(tsp_ref);
p_out1 = funcs_leak.pFn(tsp1);
p_out2 = funcs_leak2.pFn(tsp2);

figure(2); clf;
plot(tsp_ref, p_outref);
hold on
plot(tsp1, p_out1);
plot(tsp2, p_out2, '--');
hold off

%% Data
% Plot hydrophone data explicitly, multiplying by r1 (incident wave travel
% distance) to scale data near the peak pressure
floorDepth = 27;
lateralSeparation = 6;
depth = metadata_leak.paramAirgun.airgunDepth;
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

%% Wall plotting
p_L_fn = @(q) (0.4) * (q(3,:) - 0.5 * q(2,:).^2 ./ q(1,:));
T_fn = @(q) (q(3,:) - 0.5 * q(2,:).^2 ./ q(1,:)) ./ (q(1,:) * 718);

figure(5); clf;
subplot(1,2,1);
plot(solution_leak_rr0.soln.x, p_L_fn(solution_leak_rr0.q(1:3,:)));
hold on
plot(solution_leak_rr1.soln.x, p_L_fn(solution_leak_rr1.q(1:3,:)));
plot(solution_leak_rr2.soln.x, p_L_fn(solution_leak_rr2.q(1:3,:)));
plot(solution_leak_rr3.soln.x, p_L_fn(solution_leak_rr3.q(1:3,:)));

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

plot(1e-3*fieldData.t, 1e6*fieldData.y, 'k.-');
xlim([0,0.3])
legend(["\eta=1.0", "\eta=0.75", "\eta=0.5", "\eta=0.25"])
title("Flow efficiency \eta > 0.5")

%% Hydrophone data
subplot(1,2,2)

floorDepth = 27;
lateralSeparation = 6;
depth = metadata_leak_rr0.paramAirgun.airgunDepth;
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

r1_96 = norm([6, depth-9]);
hold on
plot(tsp_rr0, p_out_rr0/1e5*r1_96);
plot(tsp_rr1, p_out_rr1/1e5*r1_96);
plot(tsp_rr2, p_out_rr2/1e5*r1_96);
plot(tsp_rr3, p_out_rr3/1e5*r1_96);
% plot(tsp_rr0, p_out_rr0/1e5*r1_96, 'r');
% plot(tsp_rr1, p_out_rr1/1e5*r1_96, 'g');
% plot(tsp_rr2, p_out_rr2/1e5*r1_96, 'b');
% plot(tsp_rr3, p_out_rr3/1e5*r1_96, 'm');

% Plot DAQ
for i = 1:4
    daqIdx = sortidx(i);
    r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
    r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
    r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
    
    dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
    plot(1e-3*timeDAQ(1:3000), ...
         dataTrack(i,1:3000), ...
         'LineWidth', 1);
    hold on
end

xlim([0,0.1])
legend(["\eta=1.0", "\eta=0.75", "\eta=0.5", "\eta=0.25"])
% title("Flow efficiency \eta > 0.5")
drawnow