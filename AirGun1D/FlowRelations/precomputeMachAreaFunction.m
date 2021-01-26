% Precompute subsonic mach number from area ratio
% Returns subsonic Mach number corresponding to a particular area ratio

function exportFunction = precomputeMachAreaFunction(gamma, N)

if nargin < 2
    % Use default dictionary size
    N = 2.5e4;
end

% Sample N area ratios in interval [1, inf)
areaRatioVector = linspace(0,1,N+1);
areaRatioVector = fliplr(1./areaRatioVector(2:end));
MVector = nan(size(areaRatioVector));

% Compute M = M(A/A*)
for i = 1:length(areaRatioVector)
    areaRatio = areaRatioVector(i);
    MVector(i) = fzero( @(M) ...
        ((gamma+1)/2)^(-(gamma+1)/2/(gamma-1)) * ...
        (1 + (gamma-1)/2 * M^2 )^ ...
        ((gamma+1)/2/(gamma-1)) ./ M - ...
        areaRatio, ...
        [1e-14,1-1e-14]);
end

%% Construct interpolant
    function outputM = interpolateM(inputA)
        % Input validation
        if inputA < min(areaRatioVector)
            % Suppress any areas and return closest limit value
            inputA = 1;
%             error("precomputeMachAreaFunction: Outside of precomputed "...
%                 + "area ratio domain. Try again with larger N.")
        end
        if inputA > max(areaRatioVector)
            % Use asymptotic formula (< 1e-15 absolute error with N = 2.5e4)
            outputM = 1 ./ inputA * ...
                ((gamma + 1)/2).^(-(gamma+1)/2/(gamma-1));
        else
            % Return interpolated (vector) value of M
            outputM = interp1(areaRatioVector, MVector, inputA);
        end
    end
exportFunction = @(inputA) interpolateM(inputA);

%% Test interpolant
if false
    xAxis = linspace(1,N,2*N);
%     yAxis = g(xAxis);
    plot(xAxis, yAxis);
end

end