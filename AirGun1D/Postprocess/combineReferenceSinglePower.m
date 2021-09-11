%% Combines data from the reference and power-law bubble pressure cases
% Run after runReferenceCase and runSinglePowerCase

%% Plot bubble average pressure
figure(301); clf;

% Prep data
temp.bS_reference = [fullState.bubbleStates];
temp.bS_partition = [fullState_partition.bubbleStates];
p_reference = [temp.bS_reference.p];
t_reference = [fullState.t];
p_partition = [temp.bS_partition.p];
t_partition = [fullState_partition.t];
clear temp;

tL = tiledlayout(1,4);
for i = 1:2
    if i == 1
        nexttile(tL, 1, [1,1]);
    else
        nexttile(tL, 2, [1,3]);
    end
    plot(1e3*t_reference, p_reference/1e6, 'k', 'LineWidth', 1.5)

    % Prep field data
    psiPa_conversion = 6894.75729;
    % Voltage baseline at 0 psi
    nominalminimum = 1e5 + 9.8*10*1000;
    begin_index = 3700*5;
    end_index = begin_index+1600*5;
    fieldData = struct();
    fieldData.t = 1e3 * ...
        (HiTestData(24).iNetTimeAxisP(begin_index:end_index) - ...
         HiTestData(24).iNetTimeAxisP(begin_index));
    fieldData.y = 1e-6* (nominalminimum + ...
        psiPa_conversion * ...
        HiTestData(24).iNetCh10Data(begin_index:end_index));
    hold on
    plot(fieldData.t, fieldData.y, ...
         'LineWidth', 1.5, ...
         'Color', [0.6588    0.1686    0.1686])
    hold off
    xlim([0, 300])
    
    hold on
    plot(1e3*t_partition, p_partition/1e6, '-', ...
         'LineWidth', 1.5, ...
         'Color', [0.9608    0.5294    0.9529])
    hold off

    ylabel('$p$ [MPa]', 'Interpreter', 'latex')
    xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca, 'FontSize', 14, ...
        'TickLabelInterpreter', 'latex', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'LineWidth', 1, ...
        'TickLength', 1.5*[0.0100, 0.0250]);
    if i == 1
        xlim([0, 10])
    else
        ylabel("")
        legend(...
            [ "Model", "Data", "Partition model" ], ...
            'Interpreter', 'latex', 'Location', 'northeast');
        set(gca,'YTickLabel',arrayfun(@(x) "", 1:8))
    end
end