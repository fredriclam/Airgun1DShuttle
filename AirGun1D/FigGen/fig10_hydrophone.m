% Requires field test data
figure(10); clf;

for j = 1:2
    subplot(3,1,j);
    % Plot hydrophone data explicitly, multiplying by r1 (incident wave travel
    % distance) to scale data near the peak pressure
    floorDepth = 27;
    lateralSeparation = 6;
    depth = metadata_reference.paramAirgun.airgunDepth;
    hydrophoneDepths = [6,12,3,9];
    [~, sortidx] = sort(hydrophoneDepths);
    tSignal = HiTestData(25).headerDAQ.SamplingInterval * ...
        ((1:HiTestData(25).headerDAQ.SampleCount)-1);
    signalsData = HiTestData(25).entriesDAQ;
    start_index = 1470;
    DAQGain = 8;
    DAQSens = 1e5/7.6; % Pa per V
    manual_tshift = 3.8;
    timeDAQ = 1e3*(tSignal(start_index:end)-tSignal(start_index)) ...
        + manual_tshift;
    clear dataTrack;
    
    for i = 1:4
        daqIdx = sortidx(i);
        r1 = norm([lateralSeparation, depth-hydrophoneDepths(daqIdx)]);
        r2 = norm([lateralSeparation, depth+hydrophoneDepths(daqIdx)]);%
        r3 = norm([lateralSeparation, 2*floorDepth-depth-hydrophoneDepths(daqIdx)]);
        
        dataTrack(i,:) = DAQGain*DAQSens*signalsData(daqIdx,start_index:end)/1e5*r1;
        if j == 1 % Zoom out
            plot(timeDAQ(1:20:60000), ...
                dataTrack(i,1:20:60000), ...
                'LineWidth', 1);
        else
            plot(timeDAQ(1:3000), ...
                dataTrack(i,1:3000), ...
                'LineWidth', 1);
        end
        hold on
    end
    if j == 1 % Zoom out
        xlim([0, 2000])
    else % Zoom in
        xlim([0, 15])
    end
    drawnow
    
    disp("Data (acoustic pressure) is shown on figure.")
    
    hold on
    if j == 1 % Zoom out
        tSample = linspace(0,2,2000);
    else % Zoom in
        tSample = linspace(0,0.1,2000);
    end
    modelHydrophoneSignal = nan(4,length(tSample));
    signalAcousticFns = {};
    for i = 1:length(hydrophoneDepths)
        daqIdx = sortidx(i);
        [~, p, dists, funcs] = agtools.plotSignal(...
            fullState, ...
            solution_reference, ...
            metadata_reference, ...
            hydrophoneDepths(daqIdx), ...
            lateralSeparation, ...
            tSample);
        modelHydrophoneSignal(i,:) = p;
        signalAcousticFns{i} = @(t) dists.r1 * funcs.pFn(t)/1e5;
        disp("Plot model signal, step " + i)
    end
    
    if j == 1 % Zoom out
        xlim([0, 2000])
        ylim([-3, 5])
        set(gca,'YTick',[-3 1 5])
    else % Zoom in
        ylim([0, 5])
set(gca,'YTick',[0 2.5 5])
        xlim([0, 15])
    end
    legend(...
        [ ...
        arrayfun(@(x) ""+x+" m, data", hydrophoneDepths(sortidx)), ...
        arrayfun(@(x) ""+x+" m, model", hydrophoneDepths(sortidx))
        ], ...
        ...'Interpreter', 'latex',
        'Location', 'best');
    ylabel('$r_1 \Delta p$ (bar $\cdot$ m)')
    xlabel('{\it{t}} (ms)', 'FontSize', 14, 'interpreter', 'tex')
    set(gca, 'FontSize', 14, ...
        'TickLabelInterpreter', 'tex', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'LineWidth', 1);
    hold off
    
    % Set custom colours for plot
    colScheme =  [0         0         0
        0.3137    0.3137    0.3137
        0.5020    0.5020    0.5020
        0.8000    0.8000    0.8000
        0.3490    0.1216    0.1216
        0.6588    0.1686    0.1686
        1.0000    0.0667    0.0667
        1.0000    0.4784    0.4784];
    ch = get(gca,'children');
    for i = 1:8
        set(ch(i), 'color', colScheme(i,:))
    end
end


% Hydrophone FFT
subplot(3,1,3);

% do fft(t, p * r1), spectrum for signals spherically corrected to 1m
[wAxis, modelfft] = ...
    agtools.plotFFT((timeDAQ-timeDAQ(1))*1e-3, dataTrack(3,:)*1e5, ...
    @(t) 1e5*signalAcousticFns{3}(t));

xlabel("{\it{f}} (Hz)");%, 'Interpreter', 'latex');
ylabel("SPL (dB re 1 \mu Pa)");%, 'Interpreter', 'latex');
xlim([1e-1, 1e3]);
ylim([100, 220])
legend({'Model (9 m depth)', 'Data (9 m depth)'}, ...
    ...'Interpreter', 'latex',
    'FontSize', 14);
set(gca, 'FontSize', 14, ...
    ...'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1);
grid on

ch = get(gca,'children');
set(ch(1), 'color', [0.6588    0.1686    0.1686]);