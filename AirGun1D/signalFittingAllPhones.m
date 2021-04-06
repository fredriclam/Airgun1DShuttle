% Dependencies: curve fit toolbox

% Assuming access to solution, metadata, HiTestData

% timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));

DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
% pressureDAQ6 = DAQGain*DAQSens*HiTestData(25).entriesDAQ(1,1375:15000);
% pressureDAQ12 = DAQGain*DAQSens*HiTestData(25).entriesDAQ(2,1375:15000);
% pressureDAQ3 = DAQGain*DAQSens*HiTestData(25).entriesDAQ(3,1375:15000);
% pressureDAQ9 = DAQGain*DAQSens*HiTestData(25).entriesDAQ(4,1375:15000);
pressureDAQ6 = DAQGain*DAQSens*HiTestData(17).entriesDAQ(1,1375:15000);
pressureDAQ12 = DAQGain*DAQSens*HiTestData(17).entriesDAQ(2,1375:15000);
pressureDAQ3 = DAQGain*DAQSens*HiTestData(17).entriesDAQ(3,1375:15000);
pressureDAQ9 = DAQGain*DAQSens*HiTestData(17).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(17).headerDAQ.SamplingInterval*(1:length(pressureDAQ6));

pressures{1} = pressureDAQ6;
pressures{2} = pressureDAQ12;
pressures{3} = pressureDAQ3;
pressures{4} = pressureDAQ9;

% Set zero from noise
for i = 1:4
    pressures{i} = pressures{i} - pressures{i}(1);
end

labels{1} = '6 m depth';
labels{2} = '12 m depth';
labels{3} = '3 m depth';
labels{4} = '9 m depth';

lateralSeparation = 8;

signalFn{1} = airgunShuttleSignatureTemp(...
    solution, metadata, 6, lateralSeparation);
signalFn{2} = airgunShuttleSignatureTemp(...
    solution, metadata, 12, lateralSeparation);
signalFn{3} = airgunShuttleSignatureTemp(...
    solution, metadata, 3, lateralSeparation);
signalFn{4} = airgunShuttleSignatureTemp(...
    solution, metadata, 9, lateralSeparation);

referenceSignalFn = airgunShuttleSignatureTemp(...
    solution, metadata, 9, lateralSeparation);

%% Pre-fit plotting
figure(241); clf
for i = 1:4
subplot(2,2,i);
plot(timeDAQ, signalFn{i}(timeDAQ));
hold on
plot(timeDAQ, pressures{i}, '.');
hold off
title(labels{i});
end

%% Fitting
% Data preprocessing
dataWindowIndex = find(timeDAQ>0.016, 1,'first');
timeToFit = timeDAQ(1:dataWindowIndex)';
pressureToFit = pressureDAQ(1:dataWindowIndex)';

% Assuming referenceSignalFn, timeDAQ, pressureDAQ defined
signalFit = fit( timeToFit, pressureToFit, ...
    @(a,c,x) a.*referenceSignalFn(x-c), ...
    'Lower', [0, -0.002], ...
    'Upper', [Inf, 0.001]... % Upper limit of a, c:
    ...% upper limit on c for cannot be just quiescent signal
);

%% Check model
figure(101); clf;
plot(timeDAQ, signalFit.a*referenceSignalFn(timeDAQ-signalFit.c));
hold on
plot(timeToFit, pressureToFit, '.');
hold off

%% Check 2: scale the data instead
figure(102); clf
plot(timeDAQ, referenceSignalFn(timeDAQ));
hold on
plot(timeDAQ-signalFit.c, (1/signalFit.a)*pressureDAQ, '.');
plot(timeToFit-signalFit.c, (1/signalFit.a)*pressureToFit, '.');
hold off

%% Check all:
figure(103); clf

for i = 1:4
    subplot(2,2,i);
    plot(timeDAQ, signalFn{i}(timeDAQ));
    hold on
    plot(timeDAQ-signalFit.c, (1/signalFit.a)*pressures{i}, '.');
    hold off
    title(labels{i});
end