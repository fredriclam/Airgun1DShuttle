% Function post processing monitor data
function monitor_postprocessing(monitorStates)
monitors = [monitorStates.data.monitor];
shuttleMonitors = [monitorStates.data.shuttleMonitor];

% Monitor all points (including predictor points in ODE45)
tAxisRaw = [monitors.t];

% Sort t to make some sense of the monitor data
[tAxis, perm] = sort(tAxisRaw);
monitors = monitors(perm);
shuttleMonitors = shuttleMonitors(perm);

%% Gas distribution plot
figure(2);
subplot(3,1,1);

% Generate fill coordinates, and fill plot
fill_botContour = [shuttleMonitors.opChamberRear_m];
fill_topContour = [shuttleMonitors.opChamberRear_m] ...
                + [shuttleMonitors.opChamberFront_m];
splitFillGraph(tAxis, 0*fill_botContour, fill_botContour, fill_topContour);

xlabel 't [s]'
ylabel 'mass [kg]'

legend({'Front of shuttle', 'Behind shuttle'})

%% Pressure and Mach number
subplot(3,1,2);
plot(tAxis, [shuttleMonitors.opChamberRear_p], 'Color', [0 0 120]/255, ...
    'LineWidth', 1)
hold on
plot(tAxis, [shuttleMonitors.opChamberFront_p], 'Color', [120 0 0]/255, ...
    'LineWidth', 1)
hold off
xlabel 't [s]'
ylabel 'p [Pa]'
legend({'Front of shuttle', 'Behind shuttle'})

subplot(3,1,3);
plot(tAxis, ...
    -[shuttleMonitors.opChamberFlow_M] .* ...
      sign([shuttleMonitors.opChamberFlow_massL2R]), ...
    'Color', [0 0 0]/255, ...
    'LineWidth', 1)
xlabel 't [s]'
ylabel 'M sign(flow)'
ylim([-1, 1])

% Figure repositioning
window_pos = get(gcf, 'position');
set(gcf, 'position', [300, 50, 600, 900]);

%% Thermodynamics
figure(3);

subplot(2,1,1);
% Generate fill coordinates, and fill plot
fill_botContour = [shuttleMonitors.opChamberRear_E];
fill_topContour = [shuttleMonitors.opChamberRear_E] ...
                + [shuttleMonitors.opChamberFront_E];
splitFillGraph(tAxis, 0*fill_botContour, fill_botContour, fill_topContour);
xlabel 't [s]'
ylabel 'E [J]'
legend({'Front of shuttle', 'Behind shuttle'})

subplot(2,1,2);
% Generate fill coordinates, and fill plot
fn_entropy = @(p,rho) p.*rho.^(-1.4);
opChamberRear_entropy = [shuttleMonitors.opChamberRear_m] .*...
                          fn_entropy([shuttleMonitors.opChamberRear_p], ...
                            [shuttleMonitors.opChamberRear_rho]);
opChamberFront_entropy = [shuttleMonitors.opChamberFront_m] .*...
                          fn_entropy([shuttleMonitors.opChamberFront_p], ...
                            [shuttleMonitors.opChamberFront_rho]);
fill_botContour = opChamberRear_entropy;
fill_topContour = opChamberRear_entropy ...
                + opChamberFront_entropy;
splitFillGraph(tAxis, 0*fill_botContour, fill_botContour, fill_topContour);
xlabel 't [s]'
ylabel 'Entropy as: m log (p \rho^{ -\gamma}) [J/K]'
legend({'Front of shuttle', 'Behind shuttle'})

% Figure repositioning
window_pos = get(gcf, 'position');
set(gcf, 'position', [1000, 50, 600, 700]);

end

% Plot utility
function splitFillGraph(xAxis, curveBottom, curveMid, curveTop)
fill([xAxis, fliplr(xAxis)], [curveMid, fliplr(curveTop)],...
    [120 0 0]/255)
hold on
fill([xAxis, fliplr(xAxis)], [curveBottom, fliplr(curveMid)],...
    [0 0 120]/255)
hold off
end