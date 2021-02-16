% Main script for orchestrating code runs of airgun model with shuttle.
% 
% Uses SBPlib to set up the 1D Euler domain for the interior of the airgun
% firing chamber. Adds the model for the port area between the airgun
% domain and the bubble. Adapted from WIP test_launch_script_v1.m script.

clear; clc;

%% TODO: Plot closed chamber evolution
% TODO: Implement chambers passing data
% TODO: Implement passing exceptions
% TODO: Implement use of post-processing in Chambers
% TODO: Post check of boundary conditions in the real solution time points
% monitor_postprocessing(monitorStates)

return

%% Mass Deploy

% Disable assertion when deploying!
% File path must be correct (check the following addpath is successful)
assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

nxVector = [20, 25, 30, 35, 40, 50, 60];
pool = gcp();
for i = 1:length(nxVector)
    nx_now = nxVector(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futures(i) = parfeval(pool, @airgunShuttleDeploy, 2, nx_now);
end

% Blocking fetch
% outputs = cell(1,length(nxVector));
% for i = 1:length(nxVector)
%     [completedIndex, val] = fetchNext(futures);
%     outputs{completedIndex} = val;
%     disp("Completed job " + completedIndex)
% end

%  airgunShuttleDeploy(nx)

%% Mass Deploy II
% For different effective port areas

assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

% portAreaRatioVector = 0.6:0.1:1.4;
portAreaRatioVector = [0.6111, 0.5, 0.4, 0.3, 1.4, 2.0];
portAreaRatioVectorShort = [0.6111, 0.5, 1.1, 2.0];

% Background process or parallel (requires O(10 GB) per thread)
pool = gcp('nocreate');
if isempty(pool)
pool = parpool(1);
end

for i = 1:length(portAreaRatioVector)
    portAreaRatio = portAreaRatioVector(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futuresMidClosed(i) = parfeval(pool, @airgunShuttleDeploy, 2, 40, true, ...
        struct('airgunPortAreaRatio', portAreaRatio));
end

%% Collect data
% Use separate variable names for smaller variables
% Manually specify usable data range
for i = 1:6
%     eval(sprintf('savedResults{%d} = futures(%d).OutputArguments', i, i));
    savedResultsMidClosed{i} = futuresMidClosed(i).OutputArguments;
%     eval(sprintf("save('session-%d', 'savedResults%d', '-v7.3')", i, i));
end

%% Clear on collection
clear futuresMidClosed

%% Part 2 Runs

for i = 1:length(portAreaRatioVectorShort)
    portAreaRatio = portAreaRatioVectorShort(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futuresMidOpen(i) = parfeval(pool, @airgunShuttleDeploy, 2, 40, true, ...
        struct('airgunPortAreaRatio', portAreaRatio, ...
        'midChamberMode', 'limit-vented'));
end

for i = 1:length(portAreaRatioVectorShort)
    portAreaRatio = portAreaRatioVectorShort(i);
    
    % Call function [solution, metadata] = airgunShuttleDeploy(nx)
    futuresUncoupled(i) = parfeval(pool, @airgunShuttleDeploy, 2, 40, false, ...
        struct('airgunPortAreaRatio', portAreaRatio));
end

%% Collect data

for i = 1:4
    savedResultsMidOpen{i} = futuresMidOpen(i).OutputArguments;
end

for i = 1:4
    savedResultsUncoupled{i} = futuresUncoupled(i).OutputArguments;
end

%% Clear on collection
clear futuresMidOpen
clear futuresUncoupled

%% Pick representative results
% solution = savedResults{length(savedResults)}{1};
% metadata = savedResults{length(savedResults)}{2};
% solution = savedResultsMidOpen{2}{1};
% metadata = savedResultsMidOpen{2}{2};
solution = savedResultsMidClosed{2}{1};
metadata = savedResultsMidClosed{2}{2};

%% Postprocess result
if ~exist('fullStateBest','var')
    [fullStateBest, caseKeyContext] = ...
        airgunShuttlePostprocess(solution, metadata);
else
    airgunShuttlePostprocess(solution, metadata, fullStateBest);
end

% pressureSignalFnBest = airgunShuttleSignature(solution,metadata);

%% Work with uncoupled model

% pool = gcp('nocreate');
% if isempty(pool)
%     pool = parpool(1);
% end
% for i = 1:length(portAreaRatioVector)
%     portAreaRatio = portAreaRatioVector(i);
%     
%     % Call function [solution, metadata] = airgunShuttleDeploy(nx)
%     futuresUncoupled(i) = parfeval(pool, @airgunShuttleDeploy, 2, 40, false, ...
%         struct('airgunPortAreaRatio', portAreaRatio));
% end

%% Collect data
% Use separate variable names for smaller variables
% Manually specify usable data range
% for i = 1:9
%     savedResultsUncoupled{i} = futuresUncoupled(i).OutputArguments;
% end
% 
% % Pick last result for uncoupled model
% solnUncoupled = savedResultsUncoupled{length(savedResultsUncoupled)}{1};
% metadataUncoupled = savedResultsUncoupled{length(savedResultsUncoupled)}{2};
% fullStateUncoupled = airgunShuttlePostprocess(solnUncoupled, metadataUncoupled);

%% Compute pressure signals at far field
tSample = linspace(metadata.tspan(1), 5, 15000);

% for i = 1:length(savedResults)
%     pressureSignalFn = airgunShuttleSignature(savedResults{i}{1}, ...
%         savedResults{i}{2});
%     pressureSignals{i} = pressureSignalFn(tSample);
% end
% 
% pressureSignalFnUncoupled = airgunShuttleSignature(solnUncoupled, ...
%     metadataUncoupled);
% pressureSignalUncoupled = pressureSignalFnUncoupled(tSample);

for i = 1:length(savedResultsMidClosed)
    pressureSignalFnMidClosed = airgunShuttleSignature(savedResultsMidClosed{i}{1}, ...
        savedResultsMidClosed{i}{2});
    pressureSignalsMidClosed{i} = pressureSignalFnMidClosed(tSample);
end

for i = 1:length(savedResultsMidOpen)
    pressureSignalFnMidOpen = airgunShuttleSignature(savedResultsMidOpen{i}{1}, ...
        savedResultsMidOpen{i}{2});
    pressureSignalsMidOpen{i} = pressureSignalFnMidOpen(tSample);
end

for i = 1:length(savedResultsUncoupled)
    pressureSignalFnUncoupled = airgunShuttleSignature(savedResultsUncoupled{i}{1}, ...
        savedResultsUncoupled{i}{2});
    pressureSignalsUncoupled{i} = pressureSignalFnUncoupled(tSample);
end

%% Select model signals
modelPressureSignalMidOpen = pressureSignalsMidOpen{1};
modelPressureSignalMidClosed = pressureSignalsMidClosed{1};
modelPressureSignalUncoupled = pressureSignalsUncoupled{1};
modelPressureSignalDoubleRate = pressureSignalsDoubleRate{2};

%% Postprocess result
% solution = savedResultsMidOpen{length(savedResultsMidOpen)}{1};
% metadata = savedResultsMidOpen{length(savedResultsMidOpen)}{2};

solution = savedResultsMidClosed{4}{1};
metadata = savedResultsMidClosed{4}{2};

if ~exist('fullStateBest','var')
    [fullStateBest, caseKeyContext] = ...
        airgunShuttlePostprocess(solution, metadata);
else
    airgunShuttlePostprocess(solution, metadata, fullStateBest);
end

%% Compare Received Signal
figure(140); clf;

plot(1e3*tSample, modelPressureSignalMidOpen, '-k', 'LineWidth', 1)
hold on
plot(1e3*tSample, modelPressureSignalMidClosed, '-b', 'LineWidth', 1)
plot(1e3*tSample, modelPressureSignalUncoupled, '-r', 'LineWidth', 1)
plot(1e3*tSample, modelPressureSignalDoubleRate, '-g', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.m');
hold off

legend({...
    'mid-open', ...
    'mid-closed', ...
    'no-shuttle 100 ms', ...
    'mid-closed $2\times$ op chamber flow rate', ...
    'Data (9m depth hydrophone)'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

%% Compare Received Signal II
figure(139); clf;

plot(1e3*tSample, pressureSignalsMidClosed{1}, '-k', 'LineWidth', 1)
hold on
plot(1e3*tSample, pressureSignalsMidClosed{2}, '-b', 'LineWidth', 1)
plot(1e3*tSample, pressureSignalsMidClosed{3}, '-r', 'LineWidth', 1)
plot(1e3*tSample, pressureSignalsMidClosed{4}, '-g', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.m');
hold off

legend({...
    'mid-closed 61.1\% coverage', ...
    'mid-closed 50\% coverage', ...
    'mid-closed 40\% coverage', ...
    'mid-closed 30\% coverage', ...
    'Data (9m depth hydrophone)'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

%% Compare Received Signal III
figure(138); clf;

plot(1e3*tSample, pressureSignalsMidClosed{1}, '-k', 'LineWidth', 1)
hold on
plot(1e3*tSample, pressureSignalsMidClosed{5}, '-b', 'LineWidth', 1)
plot(1e3*tSample, pressureSignalsMidClosed{6}, '-r', 'LineWidth', 1)
% plot(1e3*tSample, pressureSignalsMidClosed{4}, '-g', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.m');
hold off

legend({...
    'mid-closed 61.1\% coverage', ...
    'mid-closed 140\% coverage', ...
    'mid-closed 200\% coverage', ...
    'Data (9m depth hydrophone)'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

%% Compare Received Signal IV
figure(138); clf;

plot(1e3*tSample, pressureSignalsDoubleRate{2}, '-k', 'LineWidth', 1)
hold on
plot(1e3*tSample, pressureSignalsDoubleRate{3}, '-r', 'LineWidth', 1)
% plot(1e3*tSample, pressureSignalsMidClosed{6}, '-r', 'LineWidth', 1)
% plot(1e3*tSample, pressureSignalsMidClosed{4}, '-g', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.m');
hold off

legend({...
    'mid-closed 61.1\% coverage, $2 \times$ op cham. flow rate, bubble 600 cui', ...
    'mid-closed 61.1\% coverage, $2 \times$ op cham. flow rate, bubble 60 cui', ...
    'Data (9m depth hydrophone)'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

%% Compare far-field pressures and bubble states
figure(141); clf;

% Select model signal
modelPressureSignalMidOpen = pressureSignalsMidOpen{1};
modelPressureSignalMidClosed = pressureSignalsMidClosed{1};
modelPressureSignalUncoupled = pressureSignalsUncoupled{1};

subplot(2,2,1);
plot(1e3*tSample, modelPressureSignalMidOpen, '-k', 'LineWidth', 1)
hold on
plot(1e3*tSample, modelPressureSignalMidClosed, '-b', 'LineWidth', 1)
plot(1e3*tSample, modelPressureSignalUncoupled, '-r', 'LineWidth', 1)
hold off
xlim([0, 500]);
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

dataDAQ = HiTestData(25).entriesDAQ(4,1375:15000);
timeDAQ = HiTestData(25).headerDAQ.SamplingInterval*(1:length(dataDAQ));
hold on
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
plot(1e3*timeDAQ, DAQGain*DAQSens*dataDAQ, '.');
hold off

legend({...
    '$T_\mathrm{L}$ mid-open', ...
    '$T_\mathrm{L}$ mid-closed', ...
    '$T_\mathrm{L}$ no-shuttle 100 ms', ...
    'Data (9m depth hydrophone)'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
grid on
grid minor

%
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

figure(142);


%% Plot pressure signals "convergence"
figure(9);
subplot(1,3,1);
for i = 1:length(savedResults)
    plot(1e3*tSample, pressureSignals{i}, '-', 'LineWidth', 1);
    hold on
end
hold off
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

legendEntries = {};
for i = 1:length(savedResults)
    legendEntries{i} = sprintf('$n_x = %d$', nxVector(i));
end
legend(legendEntries, 'Interpreter', 'latex')

subplot(1,3,2);
for i = 1:length(savedResults)-1
    plot(1e3*tSample, pressureSignals{i}-pressureSignals{end}, '-', 'LineWidth', 1);
    hold on
end
hold off
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$(\Delta p)_i - (\Delta p)_\mathrm{ref}$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend(legendEntries{1:end-1}, 'Interpreter', 'latex')

subplot(1,3,3);
for i = 1:length(savedResults)-1
    errorToRef(i) = norm(pressureSignals{i}-pressureSignals{end}, 'inf');
end
plot(nxVector(1:length(savedResults)-1), errorToRef, 'k.-', 'MarkerSize', 18);
xlabel('$n_x$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\max |(\Delta p)_i - (\Delta p)_\mathrm{ref}|$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')












return


%% [NO SHIP] - Temperature figure
% Figure not shipped with paper (too much dependence on thermocouple
% instrument response time).
% Difficulty in measuring temperature at a good time-accuracy while
% being strong enough to survive shocks within the firing chamber;
% and calibration charts for the thermocouple are typically provided for
% highly subsonic flow conditions
% May suggest in the future simulations of an incoming planar shock to the
% end of the firing chamber to find the best thermocouple configuration.

figure(901); clf;

% Convert experimental data from deg F to K
F_to_K = @(F) 5/9*(F-32) + 273.15;
T_exper_K = F_to_K(HiTestData(24).iNetCh4Data);
% Manual search index input
begin_index = 3700;
end_index = begin_index+600;

% Plot thermocouple data [K] vs. t [ms]
plot(1000*(HiTestData(24).iNetTimeAxisT(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisT(begin_index)), ...
    T_exper_K(begin_index:end_index), '.'); % [K] vs [s]
hold on
new.eDS = [fullStateBest.eulerDomainStates];
new.TAll = [new.eDS.T];
new.T_L = new.TAll(1,:);
plot(1000*[fullStateBest.t], ...
    new.T_L, ...
    'k-', 'LineWidth', 1);

% newUncoupled.eDS = [fullStateUncoupled.eulerDomainStates];
% newUncoupled.TAll = [newUncoupled.eDS.T];
% newUncoupled.T_L = [newUncoupled.TAll(1,:)];
% plot(1000*[fullStateUncoupled.t], ...
%     newUncoupled.T_L, ...
%     'r-', 'LineWidth', 1);
hold off

legend({'Data (TC 005 at end)', '$T_\mathrm{L}$', '$T_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13)
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$T$ [K]', ...
    'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

%% Figure: Plot firing chamber closed-end (x = x_L) pressure
figure(10); clf;

% Prep PDE data
new.pAll = [new.eDS.p];
new.p_L = new.pAll(1,:);
newUncoupled.pAll = [newUncoupled.eDS.p];
newUncoupled.p_L = newUncoupled.pAll(1,:);
% Prep experimental data
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % Voltage baseline at 0 psi
% Manual data index input
begin_index = 3700*5;
end_index = begin_index+1600*5;

% Plot pressure [MPa] vs. t [ms]
plot(1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisP(begin_index)), ...
    1e-6* (nominalminimum + ...
    psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index))...
    , 'b-', 'LineWidth', 1);
hold on
plot(1000*[fullStateBest.t], 1e-6*new.p_L, 'k-', 'LineWidth', 1);
plot(1000*[fullStateUncoupled.t], 1e-6*newUncoupled.p_L, 'r-', 'LineWidth', 1);

hold off

xlim([0, 500]);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Data', 'Model $p_\mathrm{L}$', 'Model $p_\mathrm{L}$ no-shuttle 100 ms'}, ...
    'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

%% Mass analysis
% analysis_solution = savedResultsUncoupled{2}{1};
% analysis_metadata = savedResultsUncoupled{2}{2};
analysis_solution = solution;
analysis_metadata = metadata;

massAirgun = sum(analysis_solution.q(1:3:end,:) ...
    * analysis_metadata.paramAirgun.airgunCrossSecAreaSqInch * 0.0254^2 ...
    * analysis_metadata.discretization.schm.h);
massBubble = analysis_solution.bubble(3,:);
figure(902); clf;
plot(analysis_solution.soln.x, massAirgun, 'r')
hold on
plot(analysis_solution.soln.x, massBubble, 'b')
plot(analysis_solution.soln.x, massAirgun + massBubble, 'k')
hold off

%% Figure: Pressure contours in x-t plot @ near time
new.t = 1e3*[fullStateBest.t];
% TODO: check metadata comes from the same source
new.x = linspace(metadata.discretization.schm.u(1), ...
    metadata.discretization.schm.u(end), size(new.pAll,1));

%% Pressure contours
figure(11); clf;

subplot(1,3,1);
plot(1e-6*new.p_L, 1000*[fullStateBest.t], 'k-', 'LineWidth', 1);
hold on
try
    psiPa_conversion = 6894.75729;
    nominalminimum = 1000*psiPa_conversion; % Voltage baseline at 0 psi
    begin_index = 3695*5;
    end_index = begin_index+1600*5;
    plot(1e-6* (nominalminimum + ...
    psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index)), ...
    1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisP(begin_index)),...
    '.', 'LineWidth', 1);
end
hold off
ylabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
xlabel('$p_\mathrm{L}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)

subplot(1,3,2);
[tGrid, xGrid] = meshgrid(new.t, new.x);
contourf(xGrid(:,1:15000), tGrid(:,1:15000), 1e-6 * new.pAll(:,1:15000), ...
    'LevelStep', 0.15);

xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
cbh = colorbar('TickLabelInterpreter', 'latex', 'FontSize', 14, ...
    'location', 'northoutside');
ylabel(cbh, '$p$ [MPa]', 'FontSize', 14, 'Interpreter', 'latex')
colormap cool
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTickLabel', {'','','','',''})
% set(gca, 'position', [0.30, 0.13, 0.45, 0.8])

% Plot tracker line for next figure
tIndexSample = 5000;
colorTrackerLine = [252, 111, 23]/255;
hold on
plot([new.x(1) new.x(end)], new.t(5000)*[1 1], ...
    'Color', colorTrackerLine, 'LineWidth', 0.5)
hold off
ylimCurrent = ylim;

% Shuttle position attachment
subplot(1,3,3);
% plot(new.shuttle_position, new.t, 'k', 'LineWidth', 1)

new.sS = [fullStateBest.shuttleStates];
new.shuttle_position = [new.sS.shuttle_position]

% Dummy lines for legend
for i = 1:6
    plot(new.shuttle_position(1), new.t(1), [caseKeyContext.colorMap{i}, '-'])
    hold on
end
for i = 1:length(caseKeyContext.caseKeySwitchIndices)-1
    % Include the right boundary element too for continuity
    plotRange = caseKeyContext.caseKeySwitchIndices(i)+1: ...
        min(length(new.t), caseKeyContext.caseKeySwitchIndices(i+1)+1);
    plot(new.shuttle_position(plotRange), new.t(plotRange), ...
        [caseKeyContext.colorMap{1+caseKeyContext.caseKeyHistory(plotRange(1))}], 'LineWidth', 1);
    hold on
end
hold off

legendLabels = {'Closed', ...
    'Subsonic', ...
    'Port choked', ...
    'Chamber choked', ...
    'Chamber choked*', ...
    'Relaxation'};
legend(legendLabels, 'Interpreter', 'latex', 'location', 'eastoutside', ...
    'FontSize', 10)

set(gca, 'YTickLabel', {'','','','',''}, ...
    'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
% Sync vertical axis
ylim([0, ylimCurrent(2)])
% set(gca, 'position', [0.76, 0.13, 0.15, 0.8])
xlabel('$\xi$ [m]', 'Interpreter', 'latex', 'FontSize', 14)

% Resize axes
subplot(1,3,3);
set(gca, 'position', [0.68, 0.13, 0.15, 0.68])
subplot(1,3,1);
ylim([0, ylimCurrent(2)])
set(gca, 'position', [0.07, 0.13, 0.16, 0.68])
subplot(1,3,2);
set(gca, 'position', [0.25, 0.13, 0.40, 0.68])

set(gcf, 'position', [-1783 399 850 483])

% Plot slice figure
figure(12); clf;
plot(new.x, 1e-6*new.pAll(:,5000), 'Color', colorTrackerLine, 'LineWidth', 1);
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 17)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
ylimCurrent = ylim;
ylim([0, ylimCurrent(2)])
% xlim([new.x(1), 0]);
xlimContours = xlim;

return

%% Figure no-shuttle Pressure contours in x-t plot @ near time
figure(111); clf;
newUncoupled.t = 1e3*[fullStateUncoupled.t];
% TODO: check metadata comes from the same source
newUncoupled.x = linspace(metadataUncoupled.discretization.schm.u(1), ...
    metadataUncoupled.discretization.schm.u(end), size(newUncoupled.pAll,1));

[tGrid, xGrid] = meshgrid(newUncoupled.t, newUncoupled.x);
contourf(xGrid(:,1:10000), tGrid(:,1:10000), 1e-6 * newUncoupled.pAll(:,1:10000), ...
    'LevelStep', 0.2);

xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
colorbar('TickLabelInterpreter', 'latex', 'FontSize', 14)
colormap cool
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

% Plot tracker line for next figure
tIndexSample = 5000;
colorTrackerLine = [0.64,0.08,0.18];
hold on
plot([newUncoupled.x(1) newUncoupled.x(end)], newUncoupled.t(5000)*[1 1], ...
    'Color', colorTrackerLine, 'LineWidth', 0.5)
hold off

figure(112); clf;
plot(newUncoupled.x, 1e-6*newUncoupled.pAll(:,5000), 'Color', colorTrackerLine, 'LineWidth', 1);
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 17)
ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
ylimCurrent = ylim;
ylim([0, ylimCurrent(2)])
xlim([newUncoupled.x(1), 0]);
return


%% Figure: Pressure contours in x-t plot @ late time
figure(13); clf;
% contourf(xGrid, tGrid, 1e-6 * new.pAll, ...
%     'LevelStep', 0.2);
contourf(xGrid(:,18000:end), tGrid(:,18000:end), 1e-6 * new.pAll(:,18000:end), ...
    'LevelStep', 0.005);

xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
colorbar('TickLabelInterpreter', 'latex', 'FontSize', 14)
colormap cool
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')

%% Figure: Pressure contours in x-t plot (global)
figure(14); clf;
contourf(xGrid, tGrid, 1e-6 * new.pAll, ...
    'LevelStep', 0.5);

xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
colorbar('TickLabelInterpreter', 'latex', 'FontSize', 14)
colormap cool
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')


%% WIP: CWT
if false
    Fs = 1/dt;
    [wt,f] = cwt(pPres,Fs);
    [m,n] = size(wt);
    subplot(4,1,[3 4]);
    h = pcolor(repmat(tInterp-r/c_inf,m,1)*1000, repmat(f,1,n), abs(wt));
    shading interp
    colormap parula
    ylim([5 220]);
    currentXLim = xlim;
    xlim([0 currentXLim(2)]);
    xlabel('Time, t-r/c_\infty (ms)');
    ylabel('Frequency (Hz)');
end