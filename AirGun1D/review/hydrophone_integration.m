% Integration of hydrophone data to give V

%%
floorDepth = 27;
lateralSeparation = 6;
depth = metadata_leak_rr0.paramAirgun.airgunDepth;
hydrophoneDepths = [6,12,3,9];
[~, sortidx] = sort(hydrophoneDepths);
tSignal = HiTestData(25).headerDAQ.SamplingInterval * ...
    ((1:HiTestData(25).headerDAQ.SampleCount)-1);
signalsData = HiTestData(25).entriesDAQ;
start_index = 1470 - 70;
DAQGain = 8;
DAQSens = 1e5/7.6; % Pa per V
manual_tshift = 3.8;
timeDAQ = 1e3*(tSignal(start_index:end)-tSignal(start_index)) ...
          + manual_tshift;
clear dataTrack;

% Plot DAQ
for i = 1:4
    daqIdx = sortidx(i);
    r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
    r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
    r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
    
    dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
end

%% 
figure(1); clf;
% Compute t shifted to zero
t = 1e-3*timeDAQ(1:maxIndex);
t = t - min(t);
empirical_V = 0*4/3*pi*0.3^3;
empirical_Vdot = 10;
% delta p at hydrophone (Pa)
maxIndex = 500;
ig = dataTrack(3,1:maxIndex) * 1e5 / r1;
ambient_mean = mean(ig(1:50));
ig = ig - ambient_mean;

empirical_Vddot_lagged = 1/(1e3 / (4*pi*r1)) * ig;
empirical_Vdot_lagged = 1/(1e3 / (4*pi*r1)) * ...
    (1e-3*(timeDAQ(2)-timeDAQ(1))) * cumsum(ig);
% Compute empirical volume up to 2 parameters: + V0 + Vdot0 * t
empirical_V_lagged = 1/(1e3 / (4*pi*r1)) * ...
    (1e-3*(timeDAQ(2)-timeDAQ(1)))^2 * cumsum(cumsum(ig));
plot(t, empirical_V_lagged);
% hold on
% plot(t, empirical_V_lagged + empirical_V +  empirical_Vdot*t);
title("Data is consistent with the following V(t) with a kernel of c_0 + c_1 * t")

%% PDE analysis
empirical_V_lagged = empirical_V_lagged - min(empirical_V_lagged) + 4/3*pi*0.3^3;
empirical_R = ((3/(4*pi))*empirical_V_lagged).^(1/3);
empirical_Rdot = empirical_Vdot_lagged ./ (4 * pi * empirical_R.^2);
empirical_Rddot = (empirical_Vddot_lagged - 8*pi*empirical_R.*empirical_Rdot.^2) ...
    ./ (4 * pi * empirical_R.^2);
% Empirical bubble pressure minus p_inf, w/o damping and 1st order compressibility term
empirical_pb = 1000 * (empirical_R .* empirical_Rddot + 1.5*(empirical_Rdot).^2);
plot(t, empirical_pb)

%% Parametric bubble pressure estimation
pm_V = @(R0, Rdot0) empirical_V_lagged - empirical_V_lagged(1) ...
    + 4/3*pi*R0^3 + 4*pi*R0^2*Rdot0*t;
pm_R = @(R0, Rdot0) ((3/(4*pi))*pm_V(R0, Rdot0)).^(1/3);
pm_Vdot = @(R0, Rdot0) empirical_Vdot_lagged - empirical_Vdot_lagged(1) + ...
    4*pi*R0^2*Rdot0;
pm_Rdot = @(R0, Rdot0) pm_Vdot(R0, Rdot0) ./ (4 * pi * pm_R(R0, Rdot0).^2);
pm_Rddot = @(R0, Rdot0) (empirical_Vddot_lagged - 8*pi*pm_R(R0, Rdot0).*pm_Rdot(R0, Rdot0).^2) ...
    ./ (4 * pi * pm_R(R0, Rdot0).^2);
pm_pb = @(R0, Rdot0) 1000 * ( ...
    pm_R(R0, Rdot0) .* pm_Rddot(R0, Rdot0) + 1.5*(pm_Rdot(R0, Rdot0)).^2);

figure(3); clf;
subplot(2,1,1)
plot(t, pm_pb(0.05, 0))
hold on
plot(t, pm_pb(0.05, 0.1))
plot(t, pm_pb(0.05, 0.5))
plot(t, pm_pb(0.05, 2.5))
for zzz = linspace(0,10.0,20)
    plot(t, pm_pb(0.05, zzz), 'k:')
end
legend(["Rdot = 0", "Rdot = 0.1", "Rdot = 0.5", "Rdot = 2.5"])

title ("Given (R=0.05, Rdot) at t = 0, p_b - p_{inf} - p_{b,0} consistent with data")


subplot(2,1,2)
plot(t, pm_pb(0.03, 0), 'LineWidth', 1.5)
hold on
plot(t, pm_pb(0.05, 0), 'LineWidth', 1.5)
plot(t, pm_pb(0.25, 0), 'LineWidth', 1.5)
plot(t, pm_pb(0.5,  0), 'LineWidth', 1.5)

for zzz = linspace(0.04,1.0,20)
    plot(t, pm_pb(zzz,  0), 'k:')
end

legend(["R = 0.03", "R = 0.05", "R = 0.25", "R = 0.5"])
title ("Given (R, Rdot=0) at t = 0, p_b - p_{inf} - p_{b,0} consistent with data")

figure(6); clf;
plot(t, pm_V(0.05, 0))
hold on
for zzz = linspace(0.05, 1.0, 20)
    plot(t, pm_V(zzz, 0), 'k:')
end

%% Euler p
euler_p = @(E, V) (1.4 - 1) * E ./ V;
p_sim = euler_p( ...
    solution_leak_rr0.bubble(4,:), ...
    solution_leak_rr0.bubble(1,:).^3 * (4*pi/3));
t_sim = solution_leak_rr0.soln.x;
figure(4); clf;
plot(t_sim, p_sim);
xlim([0, 0.016])
title("Coupled simulation p_b")

%% TEMP

figure(2); clf;
V_rr0 = solution_leak_rr0.bubble(1,:).^3 * 4/3*pi;
plot(solution_leak_rr0.soln.x, V_rr0);