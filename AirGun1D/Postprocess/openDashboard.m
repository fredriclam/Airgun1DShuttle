% Open dashboard-style figures for solution postprocessing
% Visualizes the key quantities of interest for run and typical
% quantities that can be compared to field data. 
% 
% Runs buildState to compute the postprocess data
%   Runs fullState to compute everything of interest (slow)
% 
% 
% Operates on `solution` and `metadata`, as well as `HiTestData` if it
% exists.

%% Compute postprocess data (signals, fft of signal, wall pressure)
if exist('HiTestData')
    postprocessdata = buildState( ...
        solution, metadata, HiTestData);
else
    postprocessdata = buildState( ...
        solution, metadata);
end

%% Plot mass flow rate
figure(1010); clf
tL = tiledlayout(2,2);
nexttile(4, [1, 1]);

plot(1e3*postprocessdata.t,postprocessdata.portMassFlow, '-', 'LineWidth', 1.5, ...
    'Color', [0.85, 0.33, 0.1]);
xlim([0, 300])
ylim([0, 600])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12.5, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', [0.0100    0.0250]);
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex')

legend({'Coupled model', 'Uncoupled model'}, ...
    'Location', 'best', 'Interpreter' ,'latex')

%% Plot signals
nexttile(1, [1,2]);
plot(1e3*postprocessdata.tSample, postprocessdata.signalModel/1e3, '-', 'LineWidth', 1.5);
hold on
if exist('HiTestData')
    plot(1e3*postprocessdata.timeDAQ, postprocessdata.signalData/1e3, '.-k', ...
        'MarkerSize', 2);
    
end
hold off
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\Delta p$ [kPa]', 'Interpreter', 'latex')
legend({'Coupled model', 'Hydrophone'}, ...
    'Location', 'best', 'Interpreter' ,'latex')

nexttile(3, [1, 1]);
plot(1e3*postprocessdata.tSample, postprocessdata.signal/1e3, '-', 'LineWidth', 1.5);
hold on
if exist('HiTestData')
    plot(1e3*postprocessdata.timeDAQ, postprocessdata.signalData/1e3, '.-k', ...
    'MarkerSize', 2);
end
hold off
xlim([0, 300])
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\Delta p$ [kPa]', 'Interpreter', 'latex')
legend({'Coupled model', 'Hydrophone'}, ...
    'Location', 'best', 'Interpreter' ,'latex')

%% Figure formatting
set(gcf, 'position', [-1707, 248, 631, 633]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12.5, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', [0.0100    0.0250]);

%% Frequency domain plot
figure(110); clf
loglog(postprocessdata.omegaVecPositive, ...
    abs(postprocessdata.modelfftPositive), '.-b', ...
    'MarkerSize', 2)
if exist('HiTestData')
    hold on
    loglog(postprocessdata.omegaVecDataPositive, ...
        abs(postprocessdata.datafftPositive), '.-k', ...
        'MarkerSize', 2)
    hold off
end
legend({'Coupled model', 'Data'})

%% Wall pressure
figure(301); clf;

% Plot pressure [MPa] vs. t [ms]
plot(1000*postprocessdata.t, 1e-6*postprocessdata.p_LGrid, 'k-', 'LineWidth', 1);
hold on
if exist('HiTestData')
    plot(1000*postprocessdata.iNetTimeAxis, ...
         1e-6*postprocessdata.iNetp_L, 'b-', 'LineWidth', 1);
end
hold off
xlim([0, 500]);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Model $p_\mathrm{L}$', 'Data'}, ...
    'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')