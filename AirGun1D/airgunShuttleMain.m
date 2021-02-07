

% Main script for orchestrating code runs of airgun model with shuttle.
% 
% Uses SBPlib to set up the 1D Euler domain for the interior of the airgun
% firing chamber. Adds the model for the port area between the airgun
% domain and the bubble. Adapted from WIP test_launch_script_v1.m script.

clear; clc;
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Simulation controls
% Simulation window [s]
% Suggested values:
% 0.100 to 0.600 s
tspan = [0; 0.200];
% Set flag for running shuttle-free model
runShuttleFreeFlag = false;

%% Set ambient parameters
nx = 20;                % Default: 100           % Number of grid points per 1 m of air gun length
% r = 10;                                           % Distance from source to receiver [m]
% c_inf = 1482;                                     % Speed of sound in water [m/s]
% rho_inf = 1000;                                   % Density of water [kg/m^3]
airgunDepth = 10;                                 % Depth of airgun [m]
bubbleInitialVolume = 600;                        % Initial volume [cui]
airgunPressure = 1000;                            % Initial pressure in airgun [psi]

%% Set firing chamber parameters
airgunVolume = 20600;                             % Volume of airgun [cui]
airgunInnerDiameter = 10.020;                     % Inner diameter of airgun [in]
airgunCrossSecArea = pi*airgunInnerDiameter^2/4;  % Firing chamber cross-sectional area [in^2]
airgunLength = airgunVolume / ...
    airgunCrossSecArea * 0.0254;                  % Firing chamber length [m]
% Function prescribing firing chamber profile. Not used in current version.
airgunFiringChamberProfile = @(x) error(...
    'Not implemented. Placeholder for firing chamber profile function.');

%% Set port parameters
airgunPortAreaRatio = 0.5;                        % Portion of lateral area covered by port [-]
airgunOuterDiameter = 11.2;                       % Outer diameter of firing chamber [in]
airgunPortLength = 2.5;                           % Length of port [in]
airgunPortArea = airgunPortAreaRatio * ...
    pi * airgunOuterDiameter * airgunPortLength;  % Effective port area [in^2]
% Shuttle parameters
shuttleBdryPenaltyStrength = 1e11;                % Linear elastic penalty term for shuttle [N/m]

%% Set operating chamber specifications
airCushionLength = 0.542*0.0254;                  % Length of closed air cushioning effect [m]
accelerationLength = (3.009-0.542)*0.0254;        % Length over which shuttle accelerates freely [m]

% Compression factor of air cushion as function of shuttle position
airgunOperatingChamberProfile = @(xi) (xi - accelerationLength < 0) * 1 ...
    + (xi - accelerationLength > 0) * ...
    (airCushionLength / (airCushionLength - (xi - accelerationLength)));

%% Housekeeping: Data packing
paramAirgun = struct(...
    'airgunPressure', airgunPressure, ...
    'airgunLength', airgunLength, ...
    'airgunCrossSecArea', airgunCrossSecArea, ...
    'airgunPortArea', airgunPortArea, ....
    'airgunDepth', airgunDepth, ...
    'airgunFiringChamberProfile', airgunFiringChamberProfile, ...
    'airgunOperatingChamberProfile', airgunOperatingChamberProfile, ...
    'bubbleInitialVolume', bubbleInitialVolume, ...
    'shuttleBdryPenaltyStrength', shuttleBdryPenaltyStrength);
% Initialize metadata struct for documenting results
metadata = struct(...
    'paramAirgun', paramAirgun, ...
    'nx', nx, ...
    'tspan', tspan);

%% Run solve for both models    
[solution, metadata, solShuttleFree] = ...
    runEulerCodeShuttleDual(nx, tspan, ...
                            paramAirgun, runShuttleFreeFlag, metadata);

%% Postprocess 1: Plot closed chamber evolution
% TODO: Implement chambers passing data
% TODO: Implement passing exceptions
% TODO: Implement use of post-processing in Chambers
% TODO: Post check of boundary conditions in the real solution time points
% monitor_postprocessing(monitorStates)

return

%% Mass Deploy

% Disable assertion when deploying!
% File path must be correct (check the following addpath is successful)
assert(strcmpi('C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))

addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

nxVector = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70];
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

%% Collect data
% Use separate variable names for smaller variables
% Manually specify usable data range
for i = 1:9
    eval(sprintf('savedResults{%d} = futures(%d).OutputArguments', i, i));
%     eval(sprintf("save('session-%d', 'savedResults%d', '-v7.3')", i, i));
end

% Save data
% save('session-data', 'savedResults', '-v7.3');
% Watch out! These variables can be huge.
% for i = 1:6
%     eval(sprintf('savedResults{%d} = futures(%d).OutputArguments', i, i));
%     eval(sprintf("save('session-%d', 'savedResults%d', '-v7.3')", i, i));
% end


return

%% Postprocess best result
if ~exist('fullStateBest','var')
fullStateBest = airgunShuttlePostprocess(savedResults{length(savedResults)}{1}, ...
    savedResults{length(savedResults)}{2});
else
    airgunShuttlePostprocess(savedResults{length(savedResults)}{1}, ...
    savedResults{length(savedResults)}{2}, fullStateBest);
end


%% Compute pressure signals at far field

solution = savedResults{1}{1};
metadata = savedResults{1}{2};

tSample = linspace(metadata.tspan(1), 0.1, 1000);
for i = 1:length(savedResults)
    pressureSignalFn = airgunShuttleSignature(savedResults{i}{1}, ...
        savedResults{i}{2});
    pressureSignals{i} = pressureSignalFn(tSample);
end



%% Plot pressure signals
figure(9);
subplot(1,3,1);
for i = 1:length(savedResults)
    plot(tSample, pressureSignals{i}, '-', 'LineWidth', 1);
    hold on
end
hold off
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

legendEntries = {};
for i = 1:length(savedResults)
    legendEntries{i} = sprintf('$n_x = %d$', nxVector(i));
end
legend(legendEntries, 'Interpreter', 'latex')

subplot(1,3,2);
for i = 1:length(savedResults)-1
    plot(tSample, pressureSignals{i}-pressureSignals{end}, '-', 'LineWidth', 1);
    hold on
end
hold off
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 14)
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


%%
F_to_K = @(F) 5/9*(F-32) + 273.15;

figure(10);
subplot(1,3,3);
T_exper_K = F_to_K(HiTestData(24).iNetCh4Data);
% Manual search index
begin_index = 3700;
end_index = begin_index+600;
% Plot thermocouple data [K] vs. t [ms]
plot(1000*(HiTestData(24).iNetTimeAxisT(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisT(begin_index)), ...
    T_exper_K(begin_index:end_index), '.'); % [K] vs [s]
xlim([0, 300]);
ylim([120, 300]);
title('Closed-end temperature data');
    
%% Wave travel time (experimental data)    
figure;
begin_index = 0;
end_index = 10000;
begin_index_P = 5*begin_index;
end_index_P = begin_index_P + 5*(end_index-begin_index);

% Useful downstream pressure sensor at Ch16
%   Outside sensor at Ch10
%   Upstream (port) presusre sensor at Ch13
plot(1000*(HiTestData(24).iNetTimeAxisP), ...
    HiTestData(24).iNetCh16Data);
hold on
plot(1000*(HiTestData(24).iNetTimeAxisT), ...
    HiTestData(24).iNetCh7Data); % [K] vs [s]
xlim([3660,3780])
grid on
grid minor
    
    %% Thermocouple responses
    figure; clf
    data_selection_index = 24;
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh1Data)); % 010
    hold on
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh4Data)); % 005
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh7Data)); % 003
    % Plot response time estimate using
    % Omega J-type 60 ft./sec Air
    % 427 deg C/38 deg C
    % for 005 (0.005" diameter)
    % N.B. I think we're using CHAL-005 K-type
    % - Estimate for velocity near this point?
    % http://www.newportus.com/ppt/IRCOCHAL.html
    plot(3720*[1, 1], ylim,'--')
    plot(3720*[1, 1]+80, ylim,'--')
    legend({'010', '005', '003', 'Begin Estimate', '005 response time scale estimate'});
    xlim(3600 + [0, 1500]); % for index 25
    % xlim(20000 + [0, 1000]); % for index 25
    % xlim(17500 + [0, 1000]); % for index 26
    xlabel 't [ms]'
    ylabel 'T [K]'

%% Best experimental plot
figure(2001); clf;
subplot(2,2,1);

cv = 718;
T2 = (q2(3:3:end,:) - ...
    0.5 * q2(2:3:end,:).^2 ./ q2(1:3:end,:)) ./ q2(1:3:end,:) / cv;

offset_from_left = 1;
% T_L2 = T2(30,:);
T_L2 = T2(offset_from_left, :);

subplot(2,2,1);
plot(tAxis, T_L2, 'k-', 'LineWidth', 1);
xlim([0, 300]);
ylim([120, 300]);
title('Model');
xlabel('Time [s]')
ylabel('T [K]')

subplot(2,2,2);
T_exper_K = F_to_K(HiTestData(24).iNetCh4Data);
begin_index = 3700;
end_index = begin_index+600;
plot(1000*(HiTestData(24).iNetTimeAxisT(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisT(begin_index)), ...
    T_exper_K(begin_index:end_index), 'b-', 'LineWidth', 1); % [K] vs [s]
xlim([0, 300]);
ylim([120, 300]);
title('Thermocouple 005 vs model');
xlabel('Time [s]')
ylabel('T [K]')
hold on
plot(tAxis, T_L2, 'k-', 'LineWidth', 1);
legend({'Data', 'Model'});

subplot(2,2,3);
plot(tAxis, p2(offset_from_left,:), 'k-', 'LineWidth', 1);
xlim([0, 300]);
title('Model');
xlabel('Time [s]')
ylabel('Pressure [Pa]')

subplot(2,2,4);
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % 1000 psi
begin_index = 3700*5;
end_index = begin_index+1600*5;
plot(1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisP(begin_index)), ...
    nominalminimum + ...
    psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index)...
    , 'b-', 'LineWidth', 1); % [K] vs [s]
xlim([0, 300]);
title('Pressure gauge (data + 1000 psi)');
xlabel('Time [s]')
ylabel('Pressure [Pa]')
hold on
plot(tAxis, p2(offset_from_left,:), 'k-', 'LineWidth', 1);
legend({'Data', 'Model'});

for i = 1:4
    subplot(2,2,i);
    grid minor
    grid on
end



















return

%% Call postprocess routines (legacy)

% fullState1 = airgunShuttlePostprocess(savedResults1{1}, savedResults1{2});
% airgunShuttlePostprocess(savedResults1{1}, savedResults1{2}, fullState1);

pressureSignalFn1 = airgunShuttleSignature(savedResults1{1}, savedResults1{2});
pressureSignalFn2 = airgunShuttleSignature(savedResults3{1}, savedResults2{2});
pressureSignalFn3 = airgunShuttleSignature(savedResults3{1}, savedResults3{2});
pressureSignalFn4 = airgunShuttleSignature(savedResults4{1}, savedResults4{2});
pressureSignalFn5 = airgunShuttleSignature(savedResults5{1}, savedResults5{2});
pressureSignalFn6 = airgunShuttleSignature(savedResults6{1}, savedResults6{2});

figure(9);
tSample = linspace(metadata.tspan(1), 0.1, 1000);
subplot(1,3,1);
plot(tSample, pressureSignalFn1(tSample), '-', 'LineWidth', 1);
hold on
plot(tSample, pressureSignalFn2(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn3(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn4(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn5(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn6(tSample), '-', 'LineWidth', 1);
hold off

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\Delta p$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend

subplot(1,3,2);
e1 = norm(pressureSignalFn1(tSample)-pressureSignalFn6(tSample),'inf');
e2 = norm(pressureSignalFn2(tSample)-pressureSignalFn6(tSample),'inf');
e3 = norm(pressureSignalFn3(tSample)-pressureSignalFn6(tSample),'inf');
e4 = norm(pressureSignalFn4(tSample)-pressureSignalFn6(tSample),'inf');
e5 = norm(pressureSignalFn5(tSample)-pressureSignalFn6(tSample),'inf');
plot(tSample, pressureSignalFn1(tSample)-pressureSignalFn6(tSample), '-', 'LineWidth', 1);
hold on
plot(tSample, pressureSignalFn2(tSample)-pressureSignalFn6(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn3(tSample)-pressureSignalFn6(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn4(tSample)-pressureSignalFn6(tSample), '-', 'LineWidth', 1);
plot(tSample, pressureSignalFn5(tSample)-pressureSignalFn6(tSample), '-', 'LineWidth', 1);
hold off

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$(\Delta p)_i - (\Delta p)_\mathrm{ref}$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')
legend

subplot(1,3,3);
plot([10, 15, 20, 25, 30], [e1, e2, e3, e4, e5], 'k.-', 'MarkerSize', 18);

xlabel('$n_x$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\max |(\Delta p)_i - (\Delta p)_\mathrm{ref}|$ [Pa]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

















%% Post processing
gamma = 1.4;

[~,solDY] = deval(sol, t); % Numerical differentiation (second output arg)
% Cut into the two halves: original model and shuttle-integrated model
DY1 = solDY(1:size(solDY,1)/2,:);
DY2 = solDY(size(solDY,1)/2+1:end,:);
[pPres, R, tInterp] = computePressure(bubble, DY1, t, rho_inf, ...
    c_inf, r, airgunDepth, size(q,1));
[pPres2, R2, tInterp2] = computePressure(bubble2, DY2, t, rho_inf,  ...
    c_inf, r, airgunDepth, size(q,1));

%% Plot 1: density
figure(1001); clf;
tAxis = 1000*sol.x;
tmax = max(tAxis);

subplot(1,2,1);
rho = q(1:3:end,:);
xAxis = linspace(0,1.2,size(rho, 1))';
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho', 'LineStyle', 'none');

subplot(1,2,2);
rho2 = q2(1:3:end,:);
xAxis = linspace(0,1.2,size(rho2, 1))';
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho2', 'LineStyle', 'none', 'levelstep', 2);

hold on;
plot(shuttle2(1,:)+1, tAxis, 'g');

for i = 1:2
    subplot(1,2,i);
    title('rho [kg / m^3]');
    ylim([0, tmax]);
    caxis([0, 160]);
    colorbar;
    xlabel('Position [m]');
    ylabel('Time [ms]');
end

windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    960, 350]);
colormap bone;

%% Plot 2: pressure
figure(1002); clf;

subplot(2,2,1);
p = (gamma-1)*( q(3:3:end,:)-1/2*q(2:3:end,:).^2./q(1:3:end,:) );
xAxis = linspace(0,1.2,size(p, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, p', 'LineStyle', 'none');

subplot(2,2,2);
p2 = (gamma-1)*( q2(3:3:end,:)-1/2*q2(2:3:end,:).^2./q2(1:3:end,:) );
xAxis = linspace(0,1.2,size(p2, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, p2', 'LineStyle', 'none', 'LevelStep', 0.1e6);

hold on;
aestheticOffSet = 1.2;
subsonicRange = monitorStates(3,:) == 1;
portChokedRange = monitorStates(3,:) == 2;
chamberChokedRange = monitorStates(3,:) == 3;
plot(monitorStates(17,subsonicRange) + aestheticOffSet, ...
     1000*monitorStates(13,subsonicRange), 'b.');
plot(monitorStates(17,portChokedRange) + aestheticOffSet, ...
     1000*monitorStates(13,portChokedRange), 'r.');
plot(monitorStates(17,chamberChokedRange) + aestheticOffSet, ...
     1000*monitorStates(13,chamberChokedRange), 'g.');
 
for i = 1:2
    subplot(2,2,i);
    title('p [Pa]');
    ylim([0, tmax]);
    caxis([0, 14e6]);
    colorbar;
    xlabel('Position [m]');
    ylabel('Time [ms]');
end

windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    800, 640]);
colormap bone;

subplot(2,2,3);
plot(monitorStates(3,subsonicRange), 1e3*monitorStates(13,subsonicRange), ...
    'b.', 'LineWidth', 1);
hold on
plot(monitorStates(3,portChokedRange), 1e3*monitorStates(13,portChokedRange), ...
    'r.', 'LineWidth', 1);
plot(monitorStates(3,chamberChokedRange), 1e3*monitorStates(13,chamberChokedRange), ...
    'g.', 'LineWidth', 1);
ylim([0, tmax]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Subsonic','Port-choked','Chamber-choked'})
set(gca, 'XTickLabel')
ylabel('Time [ms]');

subplot(2,2,4);
plot(shuttle2(1,:), tAxis, 'k', 'LineWidth', 1);
xlabel('Shuttle position [m]');
ylabel('Time [ms]');


%% Boundary velocity record and shuttle dynamics
% Compute first time the shuttle closes
figure(1003); clf;

u2 = q2(2:3:end,:) ./ q2(1:3:end,:);
u_R2 = u2(end,:);

subplot(2,1,1);
plot(shuttle2(1,:), tAxis, 'LineWidth', 1);
title('Shuttle pos [m]')
xlabel('xi [m]')
ylabel('t [ms]')
hold on
plot(0.0762 *[1,1], [0, tAxis(end)], '-');
plot(accelerationLength*[1,1], [0, tAxis(end)], '--');

subplot(2,1,2);
plot(shuttle2(2,:), tAxis, 'LineWidth', 1);
title('Shuttle velocity [m/s]')
hold on
plot(u_R2, tAxis, 'k', 'LineWidth', 1);
xlabel('xi-dot or u [m/s]')
ylabel('t [ms]')
legend({'Shuttle','u_R'});

%% Bubble characteristics
figure(1004); clf;
figPos = get(gcf,'Position');
set(gcf,'Position',[figPos(1) 0 600 900]);

% bubble volume
V = 4/3*pi*R.^3;
subplot(3,1,1);
plot(t*1000,R, 'LineWidth', 1);
xlabel('t [ms]');
ylabel('Radius [m]');
hold on
plot(t*1000,R2,'k','LineWidth', 1);

subplot(3,1,2);
plot(t*1000,V, 'LineWidth', 1);
xlabel('t [ms]');
ylabel('Volume [m^3]');

hold on;
V2 = 4/3*pi*R2.^3;
plot(t*1000, V2, 'k', 'LineWidth', 1);
legend({'Instant open 10 ms','Shuttle'},'location','best')

% acoustic pressure
subplot(3,1,3);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r, 'LineWidth', 1);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
hold on;
plot((tInterp2-r/c_inf)*1000, pPres2*1e-5*r, 'k', 'LineWidth', 1);


%% Down-bore temperature msmt comparison
figure(1005); clf;
cv = 718;
T = (q(3:3:end,:) - ...
    0.5 * q(2:3:end,:).^2 ./ q(1:3:end,:)) ./ q(1:3:end,:) / cv;
T2 = (q2(3:3:end,:) - ...
    0.5 * q2(2:3:end,:).^2 ./ q2(1:3:end,:)) ./ q2(1:3:end,:) / cv;

T_L = T(30,:);
T_L2 = T2(30,:);

subplot(1,3,1);
plot(tAxis, T_L);

subplot(1,3,2);
plot(tAxis, T_L2);

for i = 1:2
    subplot(1,3,i);
    title('Closed-end temperature');
    xlabel('Time [ms]');
    ylabel('T [K]');
    xlim([0, 300]);
    ylim([120, 300]);
end

figPos = get(gcf, 'position');
set(gcf, 'position', [figPos(1:2), 1300, 420])

%% Attempt to load data
% Should not work on remote devices: need undisclosed data
data_loaded = true;
try
    load ../../LinuxShare/HiTest_Data/HiTestData_v1.mat;
catch
    disp('Failed to locate data. Ignoring exp. data plot')
    data_loaded = false;
end

F_to_K = @(F) 5/9*(F-32) + 273.15;

%% Continue down-bore temperature comparison

if data_loaded
    subplot(1,3,3);
    T_exper_K = F_to_K(HiTestData(24).iNetCh4Data);
    begin_index = 3700;
    end_index = begin_index+600;
    plot(1000*(HiTestData(24).iNetTimeAxisT(begin_index:end_index) - ...
        HiTestData(24).iNetTimeAxisT(begin_index)), ...
        T_exper_K(begin_index:end_index)); % [K] vs [s]
    xlim([0, 300]);
    ylim([120, 300]);
    title('Closed-end temperature data');
    
    %% Wave travel time (experimental data)    
    figure;
    begin_index = 0;
    end_index = 10000;
    begin_index_P = 5*begin_index;
    end_index_P = begin_index_P + 5*(end_index-begin_index);
    % Downstream pressure at Ch16
    % (Ch10: outside; Ch13: inside, upstream port saturated)
    plot(1000*(HiTestData(24).iNetTimeAxisP), ...
        HiTestData(24).iNetCh10Data);
    hold on
    plot(1000*(HiTestData(24).iNetTimeAxisT), ...
        HiTestData(24).iNetCh7Data); % [K] vs [s]
    xlim([3660,3780])
    grid on
    grid minor
    legend({'P signal at port', 'T signal at end (003)'})
    title 'One-way wave travel approx. 20 ms'
    
    %% Thermocouple responses
    figure; clf
    data_selection_index = 24;
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh1Data)); % 010
    hold on
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh4Data)); % 005
    plot(1000*(HiTestData(data_selection_index).iNetTimeAxisT), ...
        F_to_K(HiTestData(data_selection_index).iNetCh7Data)); % 003
    % Plot response time estimate using
    % Omega J-type 60 ft./sec Air
    % 427 deg C/38 deg C
    % for 005 (0.005" diameter)
    % N.B. I think we're using CHAL-005 K-type
    % - Estimate for velocity near this point?
    % http://www.newportus.com/ppt/IRCOCHAL.html
    plot(3720*[1, 1], ylim,'--')
    plot(3720*[1, 1]+80, ylim,'--')
    legend({'010', '005', '003', 'Begin Estimate', '005 response time scale estimate'});
    xlim(3600 + [0, 1500]); % for index 25
    % xlim(20000 + [0, 1000]); % for index 25
    % xlim(17500 + [0, 1000]); % for index 26
    xlabel 't [ms]'
    ylabel 'T [K]'
end

%% Cleanup
clear HiTestData;

%% Best experimental plot
figure(2001); clf;
subplot(2,2,1);

cv = 718;
T2 = (q2(3:3:end,:) - ...
    0.5 * q2(2:3:end,:).^2 ./ q2(1:3:end,:)) ./ q2(1:3:end,:) / cv;

offset_from_left = 1;
% T_L2 = T2(30,:);
T_L2 = T2(offset_from_left, :);

subplot(2,2,1);
plot(tAxis, T_L2, 'k-', 'LineWidth', 1);
xlim([0, 300]);
ylim([120, 300]);
title('Model');
xlabel('Time [s]')
ylabel('T [K]')

subplot(2,2,2);
T_exper_K = F_to_K(HiTestData(24).iNetCh4Data);
begin_index = 3700;
end_index = begin_index+600;
plot(1000*(HiTestData(24).iNetTimeAxisT(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisT(begin_index)), ...
    T_exper_K(begin_index:end_index), 'b-', 'LineWidth', 1); % [K] vs [s]
xlim([0, 300]);
ylim([120, 300]);
title('Thermocouple 005 vs model');
xlabel('Time [s]')
ylabel('T [K]')
hold on
plot(tAxis, T_L2, 'k-', 'LineWidth', 1);
legend({'Data', 'Model'});

subplot(2,2,3);
plot(tAxis, p2(offset_from_left,:), 'k-', 'LineWidth', 1);
xlim([0, 300]);
title('Model');
xlabel('Time [s]')
ylabel('Pressure [Pa]')

subplot(2,2,4);
psiPa_conversion = 6894.75729;
nominalminimum = 1000*psiPa_conversion; % 1000 psi
begin_index = 3700*5;
end_index = begin_index+1600*5;
plot(1000*(HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
    HiTestData(24).iNetTimeAxisP(begin_index)), ...
    nominalminimum + ...
    psiPa_conversion*HiTestData(24).iNetCh16Data(begin_index:end_index)...
    , 'b-', 'LineWidth', 1); % [K] vs [s]
xlim([0, 300]);
title('Pressure gauge (data + 1000 psi)');
xlabel('Time [s]')
ylabel('Pressure [Pa]')
hold on
plot(tAxis, p2(offset_from_left,:), 'k-', 'LineWidth', 1);
legend({'Data', 'Model'});

for i = 1:4
    subplot(2,2,i);
    grid minor
    grid on
end

%% Phase plot
if false
    figure(1006); clf;
    
    subplot(1,2,1);
    subsonicStates =       monitorStates(:,monitorStates(3,:)==1);
    chamberLimitedStates = monitorStates(:,monitorStates(3,:)==3);
    portLimitedStates =    monitorStates(:,monitorStates(3,:)==2);
    
    plot(subsonicStates(8,:), subsonicStates(1,:), 'b.'); hold on
    plot(chamberLimitedStates(8,:), chamberLimitedStates(1,:), 'g.');
    plot(portLimitedStates(8,:), portLimitedStates(1,:), 'r.');
    plot(monitorStates(8,:),monitorStates(1,:),'-');
    ylim([0, 2]);
    
    xlabel('$A_\mathrm{p} / A_\mathrm{cs}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$M_\mathrm{a}$', 'Interpreter', 'latex', 'FontSize', 14)
    yL = ylim;
    ylim([0, yL(2)]);
    set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
    title ('Chamber exit Mach number', ...
        'Interpreter', 'latex', 'FontSize', 12)
    
    subplot(1,2,2);
    plot(subsonicStates(8,:), subsonicStates(11,:), 'b.'); hold on
    plot(chamberLimitedStates(8,:), chamberLimitedStates(11,:), 'g.');
    plot(portLimitedStates(8,:), portLimitedStates(11,:), 'r.');
    plot(monitorStates(8,:),monitorStates(11,:),'-');
    
    xlabel('$A_\mathrm{p} / A_\mathrm{cs}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$p^* / p_\mathrm{bubble} $', 'Interpreter', 'latex', 'FontSize', 14)
    yL = ylim;
    ylim([0, yL(2)]);
    set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
    set(gca,'YScale','log')
    ylim([1e-1, 1e3])
    title ('Sonic pressure (based on $p_0$) to bubble pressure', ...
        'Interpreter', 'latex', 'FontSize', 12)
end

%% CWT
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

%% Temperature tracking
if false
    c_v = 718;
    T =((q(3:3:end,:) - 0.5 * rho .* u.^2)/c_v);
    T2 =((q2(3:3:end,:) - 0.5 * q2(1:3:end,:) .* u2.^2)/c_v);
end

%% Pressure monitoring
if false
    figure(1007); clf;
    
    yMax = 10*airgunPressure * 6894.76;
    
    % Plot areas with positive velocity
    area(monitorStates(13,:), yMax*(monitorStates(6,:) > 0), ...
        'FaceColor', [0 0.4 0], 'FaceAlpha', 0.3, 'LineStyle', 'none')
    hold on;
    plot(monitorStates(13,:), monitorStates(10,:))
    plot(monitorStates(13,:), monitorStates(15,:))
    plot(monitorStates(13,:), monitorStates(16,:))
    legend({'velocity ->', 'pFiring ->', 'pRear ->', 'pFront <-'})
    xlabel 't [s]'
    ylabel('Pressure [Pa]')
    ylim([0 yMax]);
    hold off
end

%% Console report
disp('Launch script finished.')

%% Save files to a new directory
folderNameRoot = "testLaunch";
index = 1;
while exist(folderNameRoot + sprintf('%02d',index),'dir')
    index = index + 1;
end
folderName = folderNameRoot + sprintf('%02d',index);
mkdir(folderName);
cd(folderName);
for i = 1001:1006
    figure(i);
    savefig(num2str(i));
end

save('runningData','monitorStates')
cd ..
disp("Figures saved to "+ folderName);

%% Define pressure pulse at a distance
function [pPres, R, tInterp] = computePressure(bubble, DY, t, rho_inf, c_inf, r, airgunDepth, qLength)
    R = bubble(1,:); % bubble radius [m]
    V = 4/3*pi*R.^3; % bubble volume [m^3]
    U = bubble(2,:); % bubble wall velocity [m/s]
    mass = bubble(3,:); % bubble mass [kg]
    
    % A = solDY(end-2,:); % acceleration
    A = DY(qLength+2,:); % Seek acceleration entry (located after q, and R)
    % warning('Please change above line 43, it`s not right: accessing {q, bubble} data. Consider exporting from solution method');
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r + 2*airgunDepth); %ghost

    dt = 1e-6;
    tInterp = min(tDir):dt:max(tDir); % Interpolate direct, ghost waves to same uniform grid
    pDirInterp = pchip(tDir, pDir, tInterp);
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp; % High-to-low acoustic impedance: phase reversal
end