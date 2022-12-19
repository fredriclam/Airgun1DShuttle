%% Plot energy, pressures in operating chamber
% This scripts writes to the following Figures:
%   figure(771) -- Operating chamber mass distribution
%   figure(772) -- Operating chamber pressures
%   figure(773) -- Operating chamber energy distribution
%   figure(774) -- Operating chamber system entropy per total gas mass
%   figure(775) -- Total system energy
%   figure(776) -- Total system energy, part II

% Requires:
%   data (first arg from runEulerCodeShuttleDual)
%   metadata (second arg from runEulerCodeShuttleDual)

%% Common tasks

tAxis = data.soln.x;
% Note on shuttle state vector:
%   vector z = [pos; vel; m_L; E_L; m_R; E_R]
rear.m = data.shuttle(3,:);
front.m = data.shuttle(5,:);
rear.E = data.shuttle(4,:);
front.E = data.shuttle(6,:);

% Compute thermo states
rear.T = data.shuttle(4,:) ./ data.shuttle(3,:) / ...
    metadata.discretization.physConst.c_v;
front.T = data.shuttle(6,:) ./ data.shuttle(5,:) / ...
    metadata.discretization.physConst.c_v;
rear.V = arrayfun(...
    @(xi) metadata.discretization.chambers.rearVolume(xi), ...
    data.shuttle(1,:));
front.V = arrayfun(...
    @(xi) metadata.discretization.chambers.frontVolume(xi), ...
    data.shuttle(1,:));
rear.rho = rear.m ./ rear.V;
front.rho = front.m ./ front.V;
rear.p = metadata.discretization.physConst.Q * rear.rho .* rear.T;
front.p = metadata.discretization.physConst.Q * front.rho .* front.T;

% Compute from shuttleEvolve
% (used to check the code, and also contains more useful informations
% evaluated at each point)
compare.rear.p = [];
compare.front.p = [];
for i = 1:size(data.shuttle,2)
    [~, p_rear, p_front, p_Mid, subsystemState] = shuttleEvolve( ...
        data.shuttle(:,i), 0, metadata.discretization.physConst, ...
        metadata.discretization.chambers);
    compare.rear.p(i) = p_rear;
    compare.rear.T(i) = subsystemState.opChamberRear_T;
    compare.rear.rho(i) = subsystemState.opChamberRear_rho;
    compare.front.rho(i) = subsystemState.opChamberFront_rho;
    compare.front.T(i) = subsystemState.opChamberFront_T;
    compare.front.p(i) = p_front;
    compare.L2R(i) = subsystemState.opChamberFlow_massL2R;
end

% Compute shuttle kinetic energy
shuttleKineticEnergy = ...
    0.5 * metadata.discretization.physConst.shuttleAssemblyMass * ...
    data.shuttle(2,:).^2;

% Global energy
fchamber.e = data.q(3:3:end,:);
fchamber.rho = data.q(1:3:end,:);
fchamber.u = data.q(2:3:end,:) ./ data.q(1:3:end,:);
fchamber.eKinetic = 0.5 * fchamber.rho .* fchamber.u.^2;
fchamber.eThermal = fchamber.e - fchamber.eKinetic;

% Integration with respect to H matrix (from SBP)
integH = @(x) ...
    sum(metadata.discretization.H(1:3:end,1:3:end) * x(:,:),1);
fchamber.AreaCS = metadata.paramAirgun.airgunCrossSecAreaSqInch * 0.0254^2;
fchamber.total.EKinetic = fchamber.AreaCS * integH(fchamber.eKinetic);
fchamber.total.EThermal = fchamber.AreaCS * integH(fchamber.eThermal);

bubble_energy = data.bubble(4,:);

% Mid energy
mid.p0 = 1e5 + 9.8*1e3 * 10; % Manual depth
mid.T0 = 15+273;
mid.rho0 = mid.p0 / (metadata.discretization.physConst.Q * mid.T0);
% Extract subsystem state from initial step
[~, ~, ~, ~, subsystemState] = shuttleEvolve( ...
    data.shuttle(:,1), 0, metadata.discretization.physConst, ...
    metadata.discretization.chambers); 
mid.T = mid.T0 * (subsystemState.midChamber_length ./ ...
            (subsystemState.midChamber_length-data.shuttle(1,:)));
mid.m = mid.rho0 * subsystemState.midChamber_area * ...
    subsystemState.midChamber_length ;
mid.E = mid.m * metadata.discretization.physConst.c_v * mid.T;

% Compute bubble energy rates (need to be integrated)
for i = 1:size(data.shuttle,2)
    fS = fullState( ... 
        metadata.discretization, data.q(:,i), tAxis(i), ...
        data.bubble(:,i), ...
        data.shuttle(:,i), false, false);
    bubbleInterface.dQdt(i) = fS.bubbleStates.dQdt;
    bubbleInterface.pdV(i) = fS.bubbleStates.pdV;
    bubbleInterface.dEin(i) = fS.bubbleStates.dEin;
    flangeInterface.dW(i) = fS.portStates.pPort * data.shuttle(2,i) ...
        * metadata.discretization.physConst.shuttle_area_left;
end

% Damping work rate
linearDampingWorkRate = 2.5e4*data.shuttle(2,:).^2;

% Integrate work rate and heat transfer across bubble interface
bubbleInterface.Q = 0;
bubbleInterface.work = 0;
bubbleInterface.Ein = 0;
flangeInterface.W = 0;
linearDampingWork = 0;
for i = 2:size(data.shuttle,2)
    bubbleInterface.Q(i) = bubbleInterface.Q(i-1) ...
        + 0.5*(bubbleInterface.dQdt(i) + bubbleInterface.dQdt(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
    bubbleInterface.work(i) = bubbleInterface.work(i-1) ...
        + 0.5*(bubbleInterface.pdV(i) + bubbleInterface.pdV(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
    bubbleInterface.Ein(i) = bubbleInterface.Ein(i-1) ...
        + 0.5*(bubbleInterface.dEin(i) + bubbleInterface.dEin(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
    flangeInterface.W(i) = flangeInterface.W(i-1) ...
        + 0.5*(flangeInterface.dW(i) + flangeInterface.dW(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
    linearDampingWork(i) = linearDampingWork(i-1) ...
        + 0.5*(linearDampingWorkRate(i) + linearDampingWorkRate(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
end

%% pdv tracking for mid chamber
for i = 1:size(data.shuttle,2)
%     fS = fullState( ... 
%         metadata.discretization, data.q(:,i), tAxis(i), ...
%         data.bubble(:,i), ...
%         data.shuttle(:,i), false, false);
%     bubbleInterface.dQdt(i) = fS.bubbleStates.dQdt;
%     bubbleInterface.pdV(i) = fS.bubbleStates.pdV;
%     bubbleInterface.dEin(i) = fS.bubbleStates.dEin;
%     flangeInterface.dW(i) = fS.portStates.pPort * data.shuttle(2,i) ...
%         * metadata.discretization.physConst.shuttle_area_left;

    [~, ~, ~, ~, subsystemState] = shuttleEvolve( ...
        data.shuttle(:,i), 0, metadata.discretization.physConst, ...
        metadata.discretization.chambers); 
    mid.V(i) = (subsystemState.midChamber_length - data.shuttle(1,i)) * ...
        subsystemState.midChamber_area;
    mid.p(i) = subsystemState.midChamber_p;
end

mid.pavg = (mid.p(1:end-1) + mid.p(2:end))/2;
mid.dVdt = mid.V(2:end) - mid.V(1:end-1);
mid.dW = mid.pavg .* mid.dVdt;
tAxis_dt = tAxis(2:end) - tAxis(1:end-1);

mid.W = 0;
for i = 2:size(data.shuttle,2)
    mid.W(i) = mid.W(i-1) + mid.dW(i-1);% * tAxis_dt(i-1);
end

%% Gas distribution plot
figure(771); clf

% Generate fill coordinates, and fill plot
fill_botContour = rear.m;
fill_topContour = rear.m + front.m;

% Nested tiled-layout display
tL = tiledlayout(1,4);
tLNest = tiledlayout(tL,1,1);
tLNest.Layout.Tile = 1;
tLNest.Layout.TileSpan = [1, 3];
% Plot on axes set 1
ax1 = axes(tLNest);
splitFillGraph(1e3*tAxis, ...
    0*fill_botContour, fill_botContour, fill_topContour);
massLimits = max(fill_topContour);
ylim([0, massLimits]);

% Plot shuttle position on axes set 2
ax2 = axes(tLNest);
nondimShuttlePos = data.shuttle(1,:) / ...
    metadata.discretization.chambers.total_travel_length;

% Plot dummy plot for legend on the superimposed axes set 2
splitFillGraph(tAxis(end), ...
    0*fill_botContour(end), fill_botContour(end), fill_topContour(end));
hold on
% Plot shuttle position trace, downsampled to avoid fuzzy matlab line
plot(1e3*tAxis(1:100:end), nondimShuttlePos(1:100:end), '-', ...
    'Color', [0.8, 0.8, 0.8], 'LineWidth', 2);
hold off
legend({'Front of shuttle', 'Behind shuttle', ...
        'Normalized shuttle position'}, ...
    'FontSize', 12, 'Interpreter', 'latex')

ax1.Box = 'off';
ax2.Box = 'off';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.YAxis.Limits = [0, 1];

ax1.XAxis.Label.String = '$t$ [ms]';
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.FontSize = 15;
ax1.YAxis.Label.String = '$m$ [kg]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.FontSize = 15;
ax2.YAxis.Label.String = '$\xi / \xi_\mathrm{max}$';
ax2.YAxis.Label.Interpreter = 'latex';
ax2.YAxis.Label.FontSize = 15;
ax2.XAxis.TickLabels = {};

ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 13;

ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax2.YAxis.MinorTick = 'on';

ax1.LineWidth = 1;
ax2.LineWidth = 1;

ax1.TickDir = 'out';
ax2.TickDir = 'out';
grid on

% Pick next tile in outer tiled-layout
nexttile(tL);
% Compute the linear-dimension profile of the gap between the operating
% chamber and the piston
xiRange = ...
   linspace(0, metadata.discretization.chambers.total_travel_length, 100);
profileLength = arrayfun(@(x) ...
    metadata.discretization.chambers.gapProfile(x), xiRange);
plot(1e3*profileLength, xiRange/max(xiRange), 'LineWidth', 1.5, ...
    'Color', [0, 0, 0]);
ax3 = gca;
ax3.YAxis.Limits = [0, 1];
ax3.XAxis.Limits = [0, 10];
ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = 13;
ax3.YAxis.MinorTick = 'on';
ax3.LineWidth = 1;
ax3.XAxis.Label.String = '$\Delta r$ [mm]';
ax3.XAxis.Label.Interpreter = 'latex';
ax3.XAxis.Label.FontSize = 15;
grid on

% Set figure window size
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), max(windowPos(2), 600) 925, 368]);
drawnow;

%% Pressure plot
figure(772); clf
plot(1e3*tAxis, front.p/1e6, 'Color', [120 0 0]/255, ...
    'LineWidth', 1)
hold on
% Rear pressure
plot(1e3*tAxis, rear.p/1e6, 'Color', [0 0 120]/255, ...
    'LineWidth', 1)
% Initial pressure (uniform)
initial.p = metadata.paramAirgun.airgunPressure * 6894.76;
plot(1e3*tAxis, initial.p*ones(size(tAxis))/1e6, 'k--', 'LineWidth', 1);
hold off

ylim([0, 12])

legend({'Front of shuttle', 'Behind shuttle', 'Initial'}, ...
    'FontSize', 12, 'Interpreter', 'latex')

ax1 = gca;
ax1.XAxis.Label.String = '$t$ [ms]';
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.FontSize = 15;
ax1.YAxis.Label.String = '$p$ [MPa]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.FontSize = 15;
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;
ax2.LineWidth = 1.;

% Set figure window size
drawnow;
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), max(windowPos(2), 600) 603, 373]);

%% Energy
figure(773); clf

tL = tiledlayout(1,1);

ax1 = axes(tL);
% Generate fill coordinates, and fill plot
fill_botContour = rear.E/1e3;
fill_topContour = (rear.E + front.E)/1e3;
splitFillGraph(1e3*tAxis, ...
    0*fill_botContour, fill_botContour, fill_topContour);
drawnow

ax2 = axes(tL);
splitFillGraph(1e3*tAxis(end), ...
    0*fill_botContour(end), fill_botContour(end), fill_topContour(end));
hold on
% Downsampled plot to prevent fuzzy-looking curve
plot(1e3*tAxis(1:100:end), shuttleKineticEnergy(1:100:end)/1e3, '-', ...
    'Color', [0.8, 0.8, 0.8], 'LineWidth', 2);
hold off

legend(...
    {'Front of shuttle', 'Behind shuttle', 'Shuttle kinetic energy'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 13)

ax1.Box = 'off';
ax2.Box = 'off';
ax2.Color = 'none';
ax2.YAxis.Limits = ax1.YAxis.Limits;

ax1.XAxis.Label.String = '$t$ [ms]';
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.FontSize = 15;
ax1.YAxis.Label.String = '$E$ [kJ]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.FontSize = 15;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax2.YAxis.Label.String = '';
ax2.XAxis.TickValues = [];
ax2.YAxis.TickValues = [];
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.LineWidth = 1;
ax1.TickDir = 'out';

% Set figure window size
drawnow;
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), max(windowPos(2), 600) 681, 359]);

%% Entropy growth

figure(774); clf
% Generate fill coordinates, and fill plot
fn_entropy = @(p,rho) metadata.discretization.physConst.c_v * ...
    log(p.*rho.^(-1.4));
rear.S = rear.m .* fn_entropy(rear.p, rear.rho);
front.S = front.m .* fn_entropy(front.p, front.rho);
total.S = rear.S + front.S;
total.m = rear.m + front.m;
reference.S = min(total.S);
total.s_per_mass = (total.S - reference.S) ./ total.m;

plot(1e3*tAxis, total.s_per_mass, 'k', 'LineWidth', 1);
% fill_botContour = rear.S;
% fill_topContour = rear.S + front.S;
% splitFillGraph(tAxis, 0*fill_botContour, fill_botContour, fill_topContour);

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ [ms]';
ax1.XAxis.Label.FontSize = 15;
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = "$\Delta S / m_\mathrm{tot}$ [J/(kg $\cdot$ K)]";
ax1.YAxis.Label.FontSize = 15;
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

% Set figure window size
drawnow;
windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), max(windowPos(2), 600) 610, 240]);

%% Plot energy
figure(775); clf;
kiloscale = 1e-3;
% Energy amounts
energies = { ...
     rear.E(1:100:end)*kiloscale, ...
     front.E(1:100:end)*kiloscale, ...
     shuttleKineticEnergy(1:100:end)*kiloscale, ...
     fchamber.total.EThermal(1:100:end)*kiloscale, ...
     fchamber.total.EKinetic(1:100:end)*kiloscale, ...
     bubble_energy(1:100:end)*kiloscale, ...
     ...bubbleInterface.Ein(1:100:end)*kiloscale, ...
     (-min(bubbleInterface.Q(1:100:end))+ ...
         bubbleInterface.Q(1:100:end))*kiloscale, ...
     (bubbleInterface.work(1:100:end))*kiloscale, ...
     ...0*mid.E(1:100:end)*kiloscale, ...
     (max(flangeInterface.W(1:100:end))-flangeInterface.W(1:100:end))...
         *kiloscale, ...
     linearDampingWork(1:100:end)*kiloscale, ...
     -mid.W(1:100:end)*kiloscale ...
    };
labels = {'Region III', 'Region IV', 'Shuttle KE', 'Region I thermal', ...
    'Region I kinetic', 'Bubble internal', 'Water thermal $\Delta Q + Q_0$', ...
    'Bubble work', ...
    'Neglected Region I work', 'Linear damping loss', 'Region II energy'};

% Sort energies by max magnitude
energymaxmagnitudes = nan(length(energies),1);
for i = 1:length(energies)
    energymaxmagnitudes(i) = max(abs(energies{i}));
end
[~, sortidx] = sort(energymaxmagnitudes, 'descend');

% Swap regions III and IV
sortidx(sortidx<3) = 3-sortidx(sortidx<3);
% Bring kinetic to first
sortidx(sortidx == 5) = [];
sortidx = [5; sortidx];

for i = 1:length(sortidx)
    energiesSorted{i} = energies{sortidx(i)};
    labelsSorted{i} = labels{sortidx(i)};
end

splitFillMulti(1e3*tAxis(1:100:end), ...
    energiesSorted);
legend(labelsSorted, 'location', 'eastoutside', 'Interpreter', 'latex', ...
    'FontSize', 13)

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ [ms]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ [kJ]';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), min(windowPos(2), 500) 895   458]);
drawnow;

% Part II (zoom in on small energies)
figure(776); clf;

splitFillMulti(1e3*tAxis(1:100:end), ...
    energiesSorted(6:end), 6-1);
legend(labelsSorted(6:end), 'location', 'eastoutside', ...
    'Interpreter', 'latex', ...
    'FontSize', 13)

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ [ms]';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ [kJ]';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

windowPos = get(gcf,'position');
set(gcf,'position', [windowPos(1), min(windowPos(2), 500) 740   458]);
drawnow;


%% Plot utility
function splitFillGraph(xAxis, curveBottom, curveMid, curveTop)
fill([xAxis, fliplr(xAxis)], [curveMid, fliplr(curveTop)],...
    [120 0 0]/255)
hold on
fill([xAxis, fliplr(xAxis)], [curveBottom, fliplr(curveMid)],...
    [0 0 120]/255)
hold off
end

function splitFillMulti(xAxis, cellsCurves, colorSkipIdx)
    if nargin < 3
        colorSkipIdx = 0;
    end
    colortriad = [204 166 1;
                  0 155 120;
                  120 147 245;
                  173 245 95;
                  123 168 74;
                  0   26  120;
                  120 26  0;
                  199 0 95;
                  255 168 74;
                  120 249 245;
                  165 49 214] / 255;
    
    curveLower = 0*cellsCurves{1};
    curveUpper = curveLower;
    for i = 1:length(cellsCurves)
        colorIdx = 1+mod(colorSkipIdx+i-1, size(colortriad,1));
        curveUpper = curveUpper + cellsCurves{i};
        fill([xAxis, fliplr(xAxis)], ...
             [curveLower, fliplr(curveUpper)], ...
             colortriad(colorIdx,:));
        curveLower = curveLower + cellsCurves{i};
        if i == 1
            hold on
        end
        if i == length(cellsCurves)
            hold off
        end
    end

end