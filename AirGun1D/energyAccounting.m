% Energy accounting from full state history

tVec = solution.soln.x;

%% Unpacking states for energy accounting
energyAcct.eDS = [fullStateBest.eulerDomainStates];
energyAcct.bS = [fullStateBest.bubbleStates];
energyAcct.sS = [fullStateBest.shuttleStates];
energyAcct.pS = [fullStateBest.portStates];

% plot(fullStateBest.t, energyAcct.eDS

%% Power entries
% Bubble boundary work, heat loss
energyAcct.Vdot = 4*pi*[energyAcct.bS.R].^2 .*[energyAcct.bS.RDot];
bubbleBoundaryWorkRate = [energyAcct.bS.p] .* energyAcct.Vdot;
bubbleHeatLossRate = 4*pi*10*4e3*([energyAcct.bS.R].^2).*([energyAcct.bS.T]-273);

bubbleBoundaryWork = 0;
bubbleHeatLoss = 0;
for i = 2:length(tVec)
bubbleBoundaryWork(i) = bubbleBoundaryWork(i-1) + ...
    mean([bubbleBoundaryWorkRate(i-1), bubbleBoundaryWorkRate(i)]) ...
    * (tVec(i) - tVec(i-1));
bubbleHeatLoss(i) = bubbleHeatLoss(i-1) + ...
    mean([bubbleHeatLossRate(i-1), bubbleHeatLossRate(i)]) * (tVec(i) - tVec(i-1));
end

%% Energy entries
chamberEnergy = sum([energyAcct.eDS.e],1) ...
    * metadata.discretization.schm.h ...
    * metadata.discretization.physConst.crossSectionalArea;
bubbleEnergy = [energyAcct.bS.E];
opChamberRearEnergy = [energyAcct.sS.opChamberRear_E];
opChamberFrontEnergy = [energyAcct.sS.opChamberFront_E];
shuttleEnergy = 0.5 ...
    * metadata.discretization.physConst.shuttleAssemblyMass...
    * [energyAcct.sS.shuttle_velocity].^2;
% Missing: middle spring energy (negligible)
% midChamberEnergy = 

sumEnergy = chamberEnergy + ...
    bubbleEnergy + ...
    opChamberRearEnergy + ...
    opChamberFrontEnergy + ...
    shuttleEnergy + ...
    bubbleBoundaryWork + ...
    bubbleHeatLoss;

%% Plot energy
figure(1234); clf
plot(tVec, sumEnergy, 'k', 'LineWidth', 2);
hold on
plot(tVec, chamberEnergy);
plot(tVec, bubbleEnergy);
plot(tVec, bubbleBoundaryWork);
plot(tVec, bubbleHeatLoss);
plot(tVec, opChamberRearEnergy);
plot(tVec, opChamberFrontEnergy);
plot(tVec, shuttleEnergy);
legend({ ...
    'Total energy', ...
    'Firing chamber energy', ...
    'Bubble energy', ...
    'Water pdv work', ...
    'Water heat gained', ...
    'op-rear energy', ...
    'op-front energy', ...
    'shuttle-kinetic energy', ...
})

%% On bubble heat transfer
figure(1235); clf
subplot(1,3,1);
plot(tVec, [energyAcct.bS.p]);
xlabel 't [s]'
ylabel 'p [Pa]'

subplot(1,3,2);
plot(tVec, [energyAcct.bS.rho]);
xlabel 't [s]'
ylabel '\rho [kg/m^3]'

subplot(1,3,3);
plot(tVec, [energyAcct.bS.T]);
hold on
plot(tVec, 273*ones(size(tVec)), '--');
xlabel 't [s]'
ylabel 'T [K]'