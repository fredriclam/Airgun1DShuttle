figure(15); clf;

%% Run pressure variation cases
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(2);
end
nx = 40;

futuresPressures_(1) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'bubbleModel', struct( ...
          'type', 'single', ...
          'M', 10, ...
          'alpha', 0.8), ...
        'airgunPressure', 600 ...
));
futuresPressures_(2) = parfeval(pool, @airgunShuttleDeploy, 2, ...
    nx, true, ...
    struct( ...
        'bubbleModel', struct( ...
          'type', 'single', ...
          'M', 10, ...
          'alpha', 0.8), ...
        'airgunPressure', 400 ...
));

%% Compute postprocess data (signals, fft of signal, wall pressure)
wait(futuresPressures_);
disp('Run complete.')
%% Memory management
solution_600 = futuresPressures_(1).OutputArguments{1};
metadata_600 = futuresPressures_(1).OutputArguments{2};
solution_400 = futuresPressures_(2).OutputArguments{1};
metadata_400 = futuresPressures_(2).OutputArguments{2};
clear futuresPressures_;

%% Get full state history
[fullState_600, ~] = ...
        airgunShuttlePostprocess( ...
        solution_600, ...
        metadata_600);
[fullState_400, ~] = ...
        airgunShuttlePostprocess( ...
        solution_400, ...
        metadata_400);

return
    
%%
figure(); clf
tL = tiledlayout(4,1);
nexttile(tL);
agtools.plotFiringChamber_exit(fullState);
hold on
agtools.plotFiringChamber_exit(fullState_600);
agtools.plotFiringChamber_exit(fullState_400);
xlabel("{\it{t}} (ms)", "Interpreter", "tex");
ylabel("\dot{m} (kg/s)", "Interpreter", "latex")

pressure_palette = [
    0, 0, 0;
    60, 71, 89;
    93, 116, 153;
    ]/255;

ch = get(gca,'children');
for i = 1:3
    set(ch(4-i), 'Color', pressure_palette(i,:))
end

chokedChamber = @(p) A_cs * (2/(2.4))^(2/.4+1) * 1.4 ...
    * p / sqrt(1.4 * 287 * metadata_reference.discretization.physConst.Tinf);
hold on
plot([fullState.t]*1e3, ones(size([fullState.t])) * chokedChamber(6.9e6), '--');
plot([fullState.t]*1e3, ones(size([fullState.t])) * chokedChamber(4.17e6), '--');
plot([fullState.t]*1e3, ones(size([fullState.t])) * chokedChamber(2.76e6), '--');
hold off

legend(["{\it{p}} = 6.90 MPa (1000 psi)", ...
        "{\it{p}} = 4.17 MPa (600 psi)", ...
        "{\it{p}} = 2.76 MPa (400 psi)"], 'Interpreter', 'tex')

%% Panel b
nexttile(tL);

tSample = linspace(0,0.3,4000);
[tSample, p, dists, funcs_ref] = agtools.plotSignal(...
    fullState, solution_reference, metadata_reference, 10, ...
    6, tSample, pressure_palette(1,:), '-');
hold on
[tSample, p, dists, funcs_600] = agtools.plotSignal(...
    fullState_600, solution_600, metadata_600, 10, ...
    6, tSample, pressure_palette(2,:), '-');
[tSample, p, dists, funcs_400] = agtools.plotSignal(...
    fullState_400, solution_400, metadata_400, 10, ...
    6, tSample, pressure_palette(3,:), '-');
hold off

xlabel("{\it{t}} (ms)", "Interpreter", "tex");

%% Panel c: shuttle position
nexttile(3);

temp.sS = [fullState.shuttleStates];
shuttle_position_ref = [temp.sS.shuttle_position];
temp.sS = [fullState_600.shuttleStates];
shuttle_position_600 = [temp.sS.shuttle_position];
temp.sS = [fullState_400.shuttleStates];
shuttle_position_400 = [temp.sS.shuttle_position];


caseInts_ref = digest_caseKey([fullState.portStates]);
caseInts_600 = digest_caseKey([fullState_600.portStates]);
caseInts_400 = digest_caseKey([fullState_400.portStates]);

t_ref = [fullState.t];
t_600 = [fullState_600.t];
t_400 = [fullState_400.t];

plot(1e3*t_ref, ...
     shuttle_position_ref, ...
     'Color', pressure_palette(1,:), ...
     'LineWidth', 1.5);
hold on
plot(1e3*t_600, ...
     shuttle_position_600, ...
     'Color', pressure_palette(2,:), ...
     'LineWidth', 1.5);
plot(1e3*t_400, ...
     shuttle_position_400, ...
     'Color', pressure_palette(3,:), ...
     'LineWidth', 1.5);
hold off

xlabel('{\it{t}} (ms)', 'Interpreter', 'tex', 'FontSize', 14)
ylabel('$\xi$ (m)', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
ylim([0,0.07])

%% Compute fft
nexttile(4);
start_index = 1470;
manual_tshift = 3.8;
tSignal = HiTestData(25).headerDAQ.SamplingInterval * ...
    ((1:HiTestData(25).headerDAQ.SampleCount)-1);
timeDAQ = 1e3*(tSignal(start_index:end)-tSignal(start_index)) ...
          + manual_tshift;
[omegaVec_ref, modelfft_ref] = agtools.plotFFT(timeDAQ*1e-3, [], ...
    @(t) funcs_ref.pFn(t)*dists.r1);
[omegaVec_600, modelfft_600] = agtools.plotFFT(timeDAQ*1e-3, [], ...
    @(t) funcs_600.pFn(t)*dists.r1);
[omegaVec_400, modelfft_400] = agtools.plotFFT(timeDAQ*1e-3, [], ...
    @(t) funcs_400.pFn(t)*dists.r1);

xlabel("{\it{t}} (ms)", "Interpreter", "tex");

%% Panel d: fft
nexttile(4);
semilogx(omegaVec_ref, 20*log10(abs(modelfft_ref)/1e-6), ...
    'Color', pressure_palette(1,:), "LineWidth", 1.5)
hold on
semilogx(omegaVec_600, 20*log10(abs(modelfft_600)/1e-6), ...
    'Color', pressure_palette(2,:), "LineWidth", 1.5)
semilogx(omegaVec_400, 20*log10(abs(modelfft_400)/1e-6), ...
    'Color', pressure_palette(3,:), "LineWidth", 1.5)
hold off

xlabel("{\it{f}} (Hz)", 'Interpreter', 'tex');
ylabel("SPL (dB re 1\muPa)", 'Interpreter', 'tex');
xlim([1e-1, 1e3]);
ylim([120, 220])
set(gca, 'FontSize', 14, ...
    'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 1, ...
    'XTick', [1e-1, 1e1, 1e3], ...
    'YTick', [120, 170, 220] ...
    );
grid on

% Converts port state string -> int
function caseInts = digest_caseKey(portStates)
for i = 1:length(portStates)
    if strcmpi('subsonic', portStates(i).caseKey)
        caseInts(i) = 0;
    elseif strcmpi('portChoked', portStates(i).caseKey)
        caseInts(i) = 1;
    else
        caseInts(i) = -1;
    end
end
end

