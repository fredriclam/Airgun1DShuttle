% Precompute subsonic mach number from pressure ratio p/p0
% Returns subsonic Mach number corresponding to a particular pressure ratio

function exportFunction = precomputeMachPressureFunction(gamma, N)

if nargin < 2
    % Use default dictionary size
    N = 2.5e4;
end

% Sample N area ratios in interval (0,1)
pressureRatioVector = linspace(0,1,N+2);
pressureRatioVector = pressureRatioVector(2:end-1);
MVector = nan(size(pressureRatioVector));

% Compute M = M(p/p0)
for i = 1:length(pressureRatioVector)
    pressureRatio = pressureRatioVector(i);
    MVector(i) = fzero( @(M) ...
        (1 + 0.5*(gamma-1) * M.^2 ) .^(-gamma/(gamma-1)) - ...
        pressureRatio, ...
        [1e-14,1-1e-14]);
end

%% Construct interpolant
    function outputM = interpolateM(inputP)
        % Input validation
        if inputP < min(pressureRatioVector)
            error("precomputeMachPressureFunction: Outside of precomputed "...
                + "pressure ratio domain. Redo precomp with larger N.")
        end
        % Return interpolated (vector) value of M
        outputM = interp1(pressureRatioVector, MVector, inputP);
    end
exportFunction = @(inputP) interpolateM(inputP);

%% Test interpolant
if false
    xAxis = linspace(1,N,2*N);
%     yAxis = g(xAxis);
    plot(xAxis, yAxis);
end

end