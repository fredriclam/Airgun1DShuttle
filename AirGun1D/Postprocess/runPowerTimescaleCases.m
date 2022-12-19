% Runs several test cases for the power-law pressure bubble model, with an
% equilibration  timescale.
% By changing @bubbleRHS, we can work with power-law pressure in the bubble
% (p ~ (r/R)^nu) where nu can go from an initial value (say, -2) to a 0 as
% t -> infty, via a timescale tau. Three ways to do this explicitly with t
% is to prescribe and exponential, a tanh, or a Gaussian function to nu. A
% time-invariant way that makes more sense is to partition energy in the
% bubble into kinetic (nonequilibrium) and internal energy, with only the
% latter being used to compute the pressure.

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
    pool = parpool(6);
end

nx = 40;

% Coupled, with power-law pressure in bubble
futures(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -7/3, ...
        'timescaleRelaxation', 0.2)));
futures(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -7/3, ...
        'timescaleRelaxation', 0.1)));
futures(3) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -2, ...
        'timescaleRelaxation', 0.2)));
futures(4) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -2, ...
        'timescaleRelaxation', 0.1)));
futures(5) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -7/3, ...
        'timescaleRelaxation', 0.3)));
futures(6) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct('bubbleModel', struct(...
        'type', 'single-power', ...
        'M', 10, ...
        'alpha', 0.8, ...
        'pressurePower', -7/3, ...
        'timescaleRelaxation', 0.05)));
%% Compute postprocess data (signals, fft of signal, wall pressure)
wait(futures);
disp('Beginning postprocessing.')
for i = 1:6
    postprocessStates(i) = buildState( ...
        futures(i).OutputArguments{1}, ...
        futures(i).OutputArguments{2}, HiTestData);
end

%% Create legend labels
legendlabels = {...
    'Hydrophone data', ...
    '$\nu = -7/3, \tau = 0.2$', ...
    '$\nu = -7/3, \tau = 0.1$', ...
    '$\nu = -2, \tau = 0.2$', ...
    '$\nu = -2, \tau = 0.1$', ...
    '$\nu = -7/3, \tau = 0.05$', ...
    '$\nu = 0$'};

%% Plot mass flow rate
figure(3010); clf
tL = tiledlayout(2,2);
nexttile(4, [1, 1]);

plot(1e3*postprocessStates(1).t,postprocessStates(1).portMassFlow, ...
    '--', 'LineWidth', 1.5);
hold on
plot(1e3*postprocessStates(2).t,postprocessStates(2).portMassFlow, ...
    '-', 'LineWidth', 1.5);
plot(1e3*postprocessStates(3).t,postprocessStates(3).portMassFlow, ...
    '--', 'LineWidth', 1.5);
plot(1e3*postprocessStates(4).t,postprocessStates(4).portMassFlow, ...
    '-', 'LineWidth', 1.5);
plot(1e3*postprocessStates(5).t,postprocessStates(5).portMassFlow, ...
    '--', 'LineWidth', 1.5);
plot(1e3*postprocessStates(6).t,postprocessStates(6).portMassFlow, ...
    '-', 'LineWidth', 1.5);
hold off
xlim([0, 300])
ylim([0, 600])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12.5, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', [0.0100    0.0250]);
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex')

legend(legendlabels(2:end), ...
       'Location', 'best', 'Interpreter' ,'latex')

%% Plot signals
nexttile(1, [1,2]);

linecolor = colormap('lines');

plot(1e3*postprocessStates(1).timeDAQ, ...
     postprocessStates(1).signalData/1e3, '.-k', ...
     'MarkerSize', 2);
hold on
plot(1e3*postprocessStates(1).tSample, ...
     postprocessStates(1).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(1,:));
plot(1e3*postprocessStates(2).tSample, ...
     postprocessStates(2).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(2,:));
plot(1e3*postprocessStates(3).tSample, ...
     postprocessStates(3).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(3,:));
plot(1e3*postprocessStates(4).tSample, ...
     postprocessStates(4).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(4,:));
plot(1e3*postprocessStates(5).tSample, ...
     postprocessStates(5).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(5,:));
plot(1e3*postprocessStates(6).tSample, ...
     postprocessStates(6).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(6,:));
hold off
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\Delta p$ [kPa]', 'Interpreter', 'latex')
legend(legendlabels, ...
       'Location', 'best', 'Interpreter' ,'latex')

nexttile(3, [1, 1]);
plot(1e3*postprocessStates(1).timeDAQ, ...
     postprocessStates(1).signalData/1e3, '.-k', ...
     'MarkerSize', 2);
hold on
plot(1e3*postprocessStates(1).tSample, ...
     postprocessStates(1).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(1,:));
plot(1e3*postprocessStates(2).tSample, ...
     postprocessStates(2).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(2,:));
plot(1e3*postprocessStates(3).tSample, ...
     postprocessStates(3).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(3,:));
plot(1e3*postprocessStates(4).tSample, ...
     postprocessStates(4).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(4,:));
plot(1e3*postprocessStates(5).tSample, ...
     postprocessStates(5).signalModel/1e3, '-', 'LineWidth', 1.5, ...
     'Color', linecolor(5,:));
plot(1e3*postprocessStates(6).tSample, ...
     postprocessStates(6).signalModel/1e3, '--', 'LineWidth', 1.5, ...
     'Color', linecolor(6,:));
hold off
xlim([0, 300])
xlabel('$t$ [ms]', 'Interpreter', 'latex')
ylabel('$\Delta p$ [kPa]', 'Interpreter', 'latex')

legend(legendlabels, ...
       'Location', 'best', 'Interpreter' ,'latex')
   
%% Figure formatting
set(gcf, 'position', [-1707, 248, 1000, 633]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12.5, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', [0.0100    0.0250]);

%% Frequency domain plot
% figure(1011); clf
% loglog(postprocessStates(1).omegaVecDataPositive, ...
%     abs(postprocessStates(1).datafftPositive), '.-k', ...
%     'MarkerSize', 2)
% hold on
% loglog(postprocessStates(1).omegaVecPositive, ...
%     abs(postprocessStates(1).modelfftPositive), '.-b', ...
%     'MarkerSize', 2)
% hold off
% legend({'Data','Coupled model'})

%% Wall pressure
figure(3012); clf;

% Plot pressure [MPa] vs. t [ms]
plot(1000*postprocessStates(1).iNetTimeAxis, ...
     1e-6*postprocessStates(1).iNetp_L, 'b-', 'LineWidth', 1);
hold on
plot(1000*postprocessStates(1).t, 1e-6*postprocessStates(1).p_LGrid, 'k-', 'LineWidth', 1);
plot(1000*postprocessStates(2).t, 1e-6*postprocessStates(2).p_LGrid, 'k-', 'LineWidth', 1);
plot(1000*postprocessStates(3).t, 1e-6*postprocessStates(3).p_LGrid, 'r-', 'LineWidth', 1);
plot(1000*postprocessStates(4).t, 1e-6*postprocessStates(4).p_LGrid, 'g-', 'LineWidth', 1);
plot(1000*postprocessStates(5).t, 1e-6*postprocessStates(5).p_LGrid, 'r-', 'LineWidth', 1);
plot(1000*postprocessStates(6).t, 1e-6*postprocessStates(6).p_LGrid, 'g-', 'LineWidth', 1);
hold off

xlim([0, 500]);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
legend(legendlabels, ...
       'Location', 'best', 'Interpreter' ,'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

%% External/bubble pressure
figure(3013); clf;

% Plot pressure [MPa] vs. t [ms]
% Trim index for data
trimIndex = find(postprocessStates(1).iNetp_B > 0.15e6,1 ,'first') -1 ;
iNetTimeAxis_p_B = postprocessStates(1).iNetTimeAxis_p_B(trimIndex:end) ...
    - postprocessStates(1).iNetTimeAxis_p_B(trimIndex);
plot(1000*iNetTimeAxis_p_B, ...
     1e-6*postprocessStates(1).iNetp_B(trimIndex:end), 'b-', 'LineWidth', 1);
hold on
plot(1000*postprocessStates(1).t, 1e-6*postprocessStates(1).bubblePressure, 'k-', 'LineWidth', 1);
plot(1000*postprocessStates(2).t, 1e-6*postprocessStates(2).bubblePressure, 'k-', 'LineWidth', 1);
plot(1000*postprocessStates(3).t, 1e-6*postprocessStates(3).bubblePressure, 'r-', 'LineWidth', 1);
plot(1000*postprocessStates(4).t, 1e-6*postprocessStates(4).bubblePressure, 'g-', 'LineWidth', 1);
plot(1000*postprocessStates(5).t, 1e-6*postprocessStates(5).bubblePressure, 'r-', 'LineWidth', 1);
plot(1000*postprocessStates(6).t, 1e-6*postprocessStates(6).bubblePressure, 'g-', 'LineWidth', 1);
hold off

xlim([0, 100]);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
legend(legendlabels, ...
       'Location', 'best', 'Interpreter' ,'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

%% Bubble radius
figure(3014); clf;

% Plot bubble volume
get_t = @(future) future.OutputArguments{1}.bubbleContinuationTime;
get_V = @(future) ...
    future.OutputArguments{1}.bubbleContinuationState(1,:);%.^3 * (4*pi/3);
for i = 1:6
    plot(1000*get_t(futures(i)), get_V(futures(i)), 'LineWidth', 1);
    hold on
end
hold off

legend(legendlabels(2:end), ...
       'Location', 'best', 'Interpreter' ,'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$V$ [m${}^3$]', 'Interpreter', 'latex', 'FontSize', 14)