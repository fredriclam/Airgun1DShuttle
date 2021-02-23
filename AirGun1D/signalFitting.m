% Dependencies: curve fit toolbox

% Assuming referenceSignalFn, timeDAQ, pressureDAQ defined (latter should
% be row vectors)

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
plot(timeDAQ+signalFit.c, (1/signalFit.a)*pressureDAQ, '.');
plot(timeToFit+signalFit.c, (1/signalFit.a)*pressureToFit, '.');
hold off