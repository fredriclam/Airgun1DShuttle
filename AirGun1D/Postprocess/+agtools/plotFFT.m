function [omegaVec, modelfft] = plotFFT(timeDAQ, signalData, modelSignalFn)

%% Crop t vectors
% Crop model tSample
% idxSampling = find(tSample >= timeDAQ(end), 1, 'first');
% if isempty(idxSampling)
%     % Crop data tSample
%     %     timeDAQ =
%     idxSampling = find(tSample >= timeDAQ(end), 1, 'first');
%     idxSampling = 2*ceil(idxSampling/2);
%     tSampleFFT(1:)
% %     error("tSample must sample at least the length of timeDAQ")
% else
%     tSampleFFT = tSample(1:idxSampling);
% end

%% Copy t vectors
tSample = timeDAQ;

% Trim to even:
if mod(length(tSample),2) ~= 0
    tSample = tSample(1:end-1);
    signalData = signalData(1:end-1);
end

range = @(x) x(end) - x(1);
dw = 1/range(tSample);
omegaVec = [0:1:(length(tSample)/2-1), ...
    -length(tSample)/2:1:-1] * dw;

% dB re 1 uPa
SPL_formula = @(p) 20*log10(p/1e-6);

modelfft = abs(fft(modelSignalFn(tSample)))/length(omegaVec);
semilogx(omegaVec, SPL_formula(modelfft), 'k', 'LineWidth', 1)
drawnow
hold on
datafft = abs(fft(signalData))/length(omegaVec);
try
    loglog(omegaVec, ...
        SPL_formula(datafft), 'b.-', ...
        'MarkerSize', 6, ...
        'LineWidth', 1); % 'Color', [123, 31, 21]/255)
    drawnow
    hold off
    legend(["Model", "Data"])
catch
    warning("Plotting failed. Returning")
end