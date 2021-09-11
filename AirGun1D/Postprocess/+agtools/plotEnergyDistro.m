% Plots energy distribution for reference case
% Uses some hard-coded parameter values valid for the reference case.

function plotEnergyDistro(solution, metadata)
% Replace legacy variable name
data = solution;

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
bubbleInterface.acousticEnergy = 0;
flangeInterface.W = 0;
linearDampingWork = 0;

[~, funcs] = agtools.sampleSignature( ...
    solution, metadata, 0, 9, 8);

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
    acousticEnergyIntegrand = ...
        (1 / (8*pi) * 1000 / 1482 )* ...
        (0.5*(funcs.VDotDotFn(tAxis(i)).^2 + ...
              funcs.VDotDotFn(tAxis(i-1)).^2));
    bubbleInterface.acousticEnergy(i) = bubbleInterface.acousticEnergy(i-1) ...
        + acousticEnergyIntegrand ...
             *(tAxis(i)-tAxis(i-1));
    flangeInterface.W(i) = flangeInterface.W(i-1) ...
        + 0.5*(flangeInterface.dW(i) + flangeInterface.dW(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
    linearDampingWork(i) = linearDampingWork(i-1) ...
        + 0.5*(linearDampingWorkRate(i) + linearDampingWorkRate(i-1)) ...
             *(tAxis(i)-tAxis(i-1));
end

K = @(R, Rdot) 4/3*pi * R.^3 * 1000 * 1.5 .* (Rdot).^2;
% Water kinetic energy from incompressible theory
waterKineticEnergyIncomp = K( ...
    solution.bubble(1,:), ...
    solution.bubble(2,:));

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
labels = {...
    'Region III energy', ...
    'Region IV energy', ...
    'Shuttle Kinetic Energy', ...
    'Firing chamber air internal energy', ...
    'Firing chamber air kinetic energy', ...
    'Bubble internal energy', ...
    'Water internal energy', ...
    'Bubble work', ...
    'Neglected Region I work', ...
    'Linear damping loss', ...
    'Region II energy'};

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
ax1.XAxis.Label.String = '$t$ (ms)';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ (kJ)';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

windowPos = get(gcf,'position');
% set(gcf,'position', [windowPos(1), min(windowPos(2), 500) 895   458]);
set(gcf,'position', [118   669   708   313]);
drawnow;

%% Part II (zoom in on small energies)
figure(776); clf;

splitFillMulti(1e3*tAxis(1:100:end), ...
    energiesSorted(6:end), 6-1);
legend(labelsSorted(6:end), 'location', 'eastoutside', ...
    'Interpreter', 'latex', ...
    'FontSize', 13)

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ (ms)';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ (kJ)';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

windowPos = get(gcf,'position');
% set(gcf,'position', [windowPos(1), min(windowPos(2), 500) 740   458]);
set(gcf,'position', [467   664   634   315]);
drawnow;

%% Long-time bubble energies
figure(777); clf;

bubbleLong.Q = bubbleInterface.Q;
bubbleLong.work = bubbleInterface.work;
bubbleLong.Ein = bubbleInterface.Ein;
bubbleLong.dQdt = bubbleInterface.dQdt;
bubbleLong.pdV = bubbleInterface.pdV;
bubbleLong.dEin = bubbleInterface.dEin;
bubbleLong.acousticEnergy = bubbleInterface.acousticEnergy;
tAxisLong = data.bubbleContinuationTime;

for i = (1+size(data.q,2)):size(data.bubbleContinuationState,2)
    fS = fullState( ... 
        metadata.discretization, ...
        data.q(:,end), ...
        data.bubbleContinuationTime(i), ...
        data.bubbleContinuationState(:,i), ...
        data.shuttle(:,end), ...
        false, false);
    bubbleLong.dQdt(i) = fS.bubbleStates.dQdt;
    bubbleLong.pdV(i) = fS.bubbleStates.pdV;
    bubbleLong.dEin(i) = fS.bubbleStates.dEin;
end

for i = (1+size(data.q,2)):size(data.bubbleContinuationState,2)
    bubbleLong.Q(i) = bubbleLong.Q(i-1) ...
        + 0.5*(bubbleLong.dQdt(i) + bubbleLong.dQdt(i-1)) ...
             *(tAxisLong(i)-tAxisLong(i-1));
    bubbleLong.work(i) = bubbleLong.work(i-1) ...
        + 0.5*(bubbleLong.pdV(i) + bubbleLong.pdV(i-1)) ...
             *(tAxisLong(i)-tAxisLong(i-1));
    bubbleLong.Ein(i) = bubbleLong.Ein(i-1) ...
        + 0.5*(bubbleLong.dEin(i) + bubbleLong.dEin(i-1)) ...
             *(tAxisLong(i)-tAxisLong(i-1));
%     bubbleLong.acousticEnergy(i) = 
    acousticEnergyIntegrand = ...
        (1 / (8*pi) * 1000 / 1482 )* ...
        (0.5*(funcs.VDotDotFn(tAxisLong(i)).^2 + ...
              funcs.VDotDotFn(tAxisLong(i-1)).^2));
    bubbleLong.acousticEnergy(i) = bubbleLong.acousticEnergy(i-1) ...
        + acousticEnergyIntegrand ...
             *(tAxisLong(i)-tAxisLong(i-1));
end

bubbleLong.E = data.bubbleContinuationState(4,:);

energiesLong = { ...
    bubbleLong.E(1:3:end)*kiloscale, ...
    bubbleLong.work(1:3:end)*kiloscale, ...
    (-min(bubbleInterface.Q(1:100:end))+ ...
      bubbleLong.Q(1:3:end))*kiloscale, ...
};

% for i = 1:3
%     plot(tAxisLong(1:3:end), energiesLong{i})
%     hold on
% end

splitFillMulti(tAxisLong(1:3:end), ...
    energiesLong(1:3),2);
xlim([0,3])

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ (s)';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ (kJ)';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

set(gcf,'position', [87   360   805   268]);

%%
figure(778); clf;
% bubbleInterface.totalWaterKinetic = ...
%     bubbleInterface.work - ...
%     bubbleInterface.acousticEnergy;

splitFillMulti(tAxisLong(1:3:end), ...
    {bubbleLong.acousticEnergy(1:3:end)*kiloscale}, ...
    12);
% legend(["Water kinetic", "Water acoustic"], 'location', 'eastoutside', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 13)

xlim([0, 3])

ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = '$t$ (s)';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = '$E$ (kJ)';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;


set(gcf,'position', [87    96   808   233]);

end

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
                  165 49 214; % END of main colours
                  0 0 0;%  
                  0 255 0;] / 255;
    
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