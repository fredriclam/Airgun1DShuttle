%% Pre-process
new.bS = [fullStateBest.bubbleStates];

%%

figure(61);
plot(new.t,[new.sS.shuttle_position])
plot(new.t,[new.bS.m])

figure(62);
plot(new.t,[new.eDS.p_R])
hold on
plot(new.t,[new.pS.pPort])
hold off

%%
figure(59); clf;
plot(new.t,[new.eDS.p_R])
sketch.p_L = [new.eDS.p];
sketch.p_L = sketch.p_L(1,:);
hold on
plot(new.t,sketch.p_L)

ylabel 'p [Pa]'; xlabel 't [ms]'
legend({'Model R', 'Model L'})

%% Change in states q_r and q_r(hat)
figure(63);
subplot(1,3,1);
plot(new.t, [1 0 0] * [new.pS.qPort], 'LineWidth', 2.5);
hold on
plot(new.t, [1 0 0] * [new.eDS.q_R], 'LineWidth', 1.5);
ylabel '\rho'
subplot(1,3,2);
plot(new.t, [0 1 0] * [new.pS.qPort], 'LineWidth', 2.5);
hold on
plot(new.t, [0 1 0] * [new.eDS.q_R], 'LineWidth', 1.5);
ylabel '\rho u'
subplot(1,3,3);
plot(new.t, [0 0 1] * [new.pS.qPort], 'LineWidth', 2.5);
hold on
plot(new.t, [0 0 1] * [new.eDS.q_R], 'LineWidth', 1.5);
ylabel 'p'

%% Change in states q_r and q_r(hat)
figure(64); clf
plot(new.t, [1 0 0] * [new.pS.qPort]/new.pS(1).qPort(1), 'LineWidth', 2.5);
hold on
plot(new.t, [0 1 0] * [new.pS.qPort]/max([0 1 0] * [new.pS.qPort]), ...
    'LineWidth', 2.5);
plot(new.t, [0 0 1] * [new.pS.qPort]/new.pS(1).qPort(3), 'LineWidth', 2.5);

%% Case?
figure(65);
plot(new.t, caseKeyContext.caseKeyHistory)

%% Pressures
figure(66);
plot(new.t, [new.eDS.p_R])
hold on
plot(new.t, [new.bS.p])

%% KEY: pressure plateau prediction
sketch.gamma = 1.4;
sketch.M0 = 0.1;
exact_pRatio = @(M) (2./((sketch.gamma-1)*M+2)).^ ...
    (2*sketch.gamma/(sketch.gamma-1));

gExponent = @(g) (g+1)/2/(g-1);
areaToSonicArea = @(g,M) ((g + 1)/2)^(-gExponent(g)) .* ...
    (1 + (g-1)/2*M.^2).^gExponent(g) ./ M;
MUpstreamFn = @(ARatio) fzero(@(M) areaToSonicArea(1.4,M) - ARatio, [1e-8, 1]);
AreaRatioToDepressurization = @(ARatios) exact_pRatio(...
    arrayfun(MUpstreamFn, ARatios));

figure(67); clf
subplot(1,2,1);
sketch.MLine = linspace(1e-3, 1 - 1e-3, 100);
plot(sketch.MLine, exact_pRatio(sketch.MLine))
xlabel 'M_R'
ylabel 'p/p_{init}'

subplot(1,2,2);
sketch.ARatios = linspace(1+1e-3, 10, 100);
plot(sketch.ARatios, AreaRatioToDepressurization(sketch.ARatios))
xlabel 'A_{cs} / A_{port}'
ylabel 'p/p_{init}'

%% Pressure plateau maximum in data
figure(68); clf;

subplot(1,2,1);
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % Voltage baseline at 0 psi
begin_index = 3695*5;
end_index = begin_index+1600*5;
sketch.tData = 1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
HiTestData(24).iNetTimeAxisP(begin_index));
sketch.pData = 1e-6* (nominalminimum + ...
psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index));
plot(sketch.tData, sketch.pData, '.', 'LineWidth', 1);
xlim([0 200]);

pInitialSensor = max(sketch.pData(1:1500));
hold on
theoreticalMaxDepressurization = exact_pRatio(1.0) * pInitialSensor;
plot([0, 200], theoreticalMaxDepressurization* [1, 1]);
xlabel 't [ms]'
ylabel 'p_L [MPa]'

subplot(1,2,2);
plot(sketch.tData, sketch.pData);
hold on
plot(sketch.tData-20, 0.5*sketch.pData+0.5*max(sketch.pData));
plot([0, 200], theoreticalMaxDepressurization*[1, 1], '--');
xlabel 't [ms]'
ylabel 'p_L [MPa]'
xlim([0 200]);

%% ~Bubble pressure sensor
% May indicate that bubble model is crucial in getting more accurate
% far-field signature.
figure(269); clf;
plot(sketch.tData(19:end) - sketch.tData(18), 1e-6* (...
psiPa_conversion*HiTestData(24).iNetCh10Data(18+begin_index:end_index)));
hold on
plot(new.t,1e-6*([new.bS.p]-new.bS(1).p))
xlabel 't [ms]'
ylabel 'P [MPa] just outside port (data) or bubble (model)'
grid minor
grid on

%% Port inside sensor
% Untrustworthy: max to min dip in one sampling period. Random artifacts in
% late time may indicate sensor failure in some way
figure(70); clf;
plot(sketch.tData, 1e-6* (nominalminimum + ...
psiPa_conversion*HiTestData(24).iNetCh13Data(begin_index:end_index)));

%% Mass rate out
% Untrustworthy: max to min dip in one sampling period. Random artifacts in
% late time may indicate sensor failure in some way
figure(71); clf;
plot(new.t, [new.pS.massFlowPort]);
xlabel 't [ms]'
ylabel '{dm/dt} [kg/s]'

%% All sensors
figure(241); clf;

% Select model signal
pressureSignalFnBubble = airgunShuttleSignature(solutionSplitBubble, metadataSplitBubble);
modelPressureSignal = pressureSignalFnBubble(tSample);

subplot(2,2,1);
plot(1e3*tSample, modelPressureSignal, '-k', 'LineWidth', 1)

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.');
hold off
%%
subplot(2,2,2);
new.bubbleStates = [solution.bubbleContinuationState];
newUncoupled.bubbleStates = [solnUncoupled.bubbleContinuationState];

plot(1e3*[solution.bubbleContinuationTime], new.bubbleStates(1,:), ...
    '-k', 'LineWidth', 1)
hold on
plot(1e3*[solnUncoupled.bubbleContinuationTime], newUncoupled.bubbleStates(1,:), '-r', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$R$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend({...
%     'Data (TC 005 at end)', ...
    '$T_\mathrm{L}$', '$T_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13, ...
    'Location', 'best')
grid on
grid minor

subplot(2,2,3);
plot(1e3*[solution.bubbleContinuationTime], new.bubbleStates(2,:), ...
    '-k', 'LineWidth', 1)
hold on
plot(1e3*[solnUncoupled.bubbleContinuationTime], newUncoupled.bubbleStates(2,:), '-r', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\dot{R}$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend({...
%     'Data (TC 005 at end)', ...
    '$T_\mathrm{L}$', '$T_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

subplot(2,2,4);
tAxis1 = 1e3*[solution.bubbleContinuationTime];
tAxis1 = tAxis1(1:end-1);
RDotDot1 = diff( new.bubbleStates(2,:) ) ./ diff (1e3*[solution.bubbleContinuationTime]);

tAxis2 = 1e3*[solnUncoupled.bubbleContinuationTime];
tAxis2 = tAxis2(1:end-1);
RDotDot2 = diff( newUncoupled.bubbleStates(2,:) ) ./ diff (1e3*[solnUncoupled.bubbleContinuationTime]);

plot(tAxis1, RDotDot1, ...
    '-k', 'LineWidth', 1)
hold on
plot(tAxis2, RDotDot2, '-r', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\ddot{R}$ [m/s${}^2$]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend({...
%     'Data (TC 005 at end)', ...
    '$T_\mathrm{L}$', '$T_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

for i = 1:4
    subplot(2,2,i);
    xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
grid on
grid minor
end