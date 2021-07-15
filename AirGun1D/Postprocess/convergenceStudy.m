%% Performs convergence study (cheap version)

%% Header
assert(strcmpi(...
    'C:\Users\Fredric\Documents\Airgun\airgun_code_li\AirGun1D', ...
    cd))
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

%% Global parameters
nx = 40;

%% Bubble comparison
% Background process or parallel (requires O(10 GB) per thread)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(4);
end

nxVector = [20,30,40,50,60,70];

for i = 1:length(nxVector)
    nx = nxVector(i);
    fConvergence(i) = parfeval(pool, @airgunShuttleDeploy, 2, ...
        nx, true, ...
        struct('bubbleModel', 'single'));
end

%% On completion: extract data
for i = 1:length(nxVector)
    dataConvergence{i} = fConvergence(i).OutputArguments{1};
    metadataConvergence{i} = fConvergence(i).OutputArguments{2};
end

%% Analysis

tMax = max(dataConvergence{1}.bubbleContinuationTime);
tFine = linspace(0,0.2,1001);
tCoarse = linspace(0.2,tMax,1001);
tSample = [tFine, tCoarse(2:end)];
dtFine = tFine(2) - tFine(1);
dtCoarse = tCoarse(2) - tCoarse(1);

% Resample to common grid
for i = 1:length(dataConvergence)
    RSampledFine{i} = interp1( ...
        dataConvergence{i}.bubbleContinuationTime, ...
        dataConvergence{i}.bubbleContinuationState(1,:),tFine,'pchip');
    RSampledCoarse{i} = interp1( ...
        dataConvergence{i}.bubbleContinuationTime, ...
        dataConvergence{i}.bubbleContinuationState(1,:),tCoarse,'pchip');
    RDotSampledFine{i} = interp1( ...
        dataConvergence{i}.bubbleContinuationTime, ...
        dataConvergence{i}.bubbleContinuationState(2,:),tFine,'pchip');
    RDotSampledCoarse{i} = interp1( ...
        dataConvergence{i}.bubbleContinuationTime, ...
        dataConvergence{i}.bubbleContinuationState(2,:),tCoarse,'pchip');
end

%% Compute error
% Trap integration
normL2 = @(f, h) sqrt(h * (sum(f.^2) - 0.5*f(1)^2 - 0.5*f(end)^2));

for i = 1:length(dataConvergence)-1
    RErrs(i) = sqrt(normL2(RSampledFine{i} - RSampledFine{end}, dtFine)^2 + ...
        normL2(RSampledCoarse{i} - RSampledCoarse{end}, dtCoarse)^2);
    RDotErrs(i) = sqrt(normL2(RDotSampledFine{i} - RDotSampledFine{end}, dtFine)^2 + ...
        normL2(RDotSampledCoarse{i} - RDotSampledCoarse{end}, dtCoarse)^2);
end

%% Plot sequential errors
figure(612); clf;
subplot(2,1,1);
loglog(nxVector(1:end-1), RErrs, 'k.-');
hold on
plot(nxVector(1:end-1), 1e-5./nxVector(1:end-1), '--')
ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = ...
    '$n_x$';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = ...
    ...'$\left( \int (R^{(k)} - R^\mathrm{(ref)})^2(t) dt \right)^{1/2}$';
    '$\|R^{(k)} - R^\mathrm{(ref)}\|_2$';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;

subplot(2,1,2);
loglog(nxVector(1:end-1), RDotErrs, 'k.-');
hold on
plot(nxVector(1:end-1), 1e-3./nxVector(1:end-1), '--')
ax1 = gca;
ax1.XAxis.Label.Interpreter = 'latex';
ax1.XAxis.Label.String = ...
    '$n_x$';
ax1.YAxis.Label.Interpreter = 'latex';
ax1.YAxis.Label.String = ...
    ...'$\left( \int (\dot{R}^{(k)} - \dot{R}^\mathrm{(ref)})^2(t) dt \right)^{1/2}$';
    '$\|\dot{R}^{(k)} - \dot{R}^\mathrm{(ref)}\|_2$';
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.XAxis.MinorTick = 'on';
ax1.YAxis.MinorTick = 'on';
ax1.LineWidth = 1.;
