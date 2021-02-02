function [qTarget, exitFlag] = enforceScalarConstraint(obj, essentialConstraint, q_R)
% q_R: 3-element vector of the existing state
% essentialConstraint: function handle for the physical constraint

% Create projection matrix to outgoing characterisics
P = [1 0 0
     0 1 0];
 
%% Multi-D fsolve with any type of eigenvector
% Legacy code for multi-dimensional fsolve
% The difficulty here is the dependence on the numerical parameters
% (tolerance) and the resulting slowness. Use of this strategy is not
% recommended.

if false
    outgoingCharBCConstraints = @(q) ...
        P ...
        * (mapq2characteristicsStatic(q) ...
        - mapq2characteristicsStatic(q_R));
    % Add scaling to equation system to give weight to essential
    % constraint
    scalingMatrix = diag([max(1,min(1e3,1/M_R)), 1e6, 1e6]);
    unscaledConstraintSystem = @(q) [
        essentialConstraint(q);
        outgoingCharBCConstraints(q);
        ];

    % Solve constraint system
    [qPort, constraintVals, exitFlag] = fsolve(...
        @(q) scalingMatrix * unscaledConstraintSystem(q), ...
        q_R, ...
        optimoptions('fsolve', ...
        'display','none', ...
        'FunctionTolerance', 1e-7, ...
        'OptimalityTolerance', 1e-10, ...
        'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 10000, ...
        'FiniteDifferenceType', 'central' ...
        )...
        );
end

%% 1-D fsolve with static eigenvectors
% Seek to compute the target state q (q-hat) as a deviation from q_R. The
% deviation should live in the null space of some vector

% Compute null space of span({outgoing characteristic vectors})
%   Deviations from q_R live in the null space of B
%   Assumes static T(q_R) in
%     T(q_R) \ q_R == T(q_R) \ qTarget
B = P / (obj.schm.T(q_R));

try
    nullVector = null(B);
catch
    warning('null failed')
    qTarget = NaN(size(q_R));
    exitFlag = -9;
    return
end
if size(nullVector,2) > 1
    error('Nullity of characteristics preservation > 1. Unexpected!')
end

% Constrain deviation from q_R to the feasible space
qConstrained = @(deviationScalar) q_R + deviationScalar*nullVector;

% 1-D solve for target state that satisfies the essential constraint
try
[deviationScalar, ~, exitFlag] = ...
    fzero(@(deviationScalar) ...
          essentialConstraint(qConstrained(deviationScalar)), ...
          0, ...
          optimset('Display','none'));
catch
    warning('fzero failed')
    deviationScalar = NaN;
    exitFlag = -9;
end

if exitFlag == -4
    try
    % Complex fallback
    [deviationScalar, ~, exitFlag] = ...
    fzero(@(deviationScalar) ...
          essentialConstraint(qConstrained(deviationScalar)), ...
          [-1e10, 1e10], ...
          optimset('Display','none'));
    catch
        % Persist exit flag
        exitFlag = -4;
    end
end
      
% Compute target state
qTarget = q_R + deviationScalar * nullVector;                  

% fsolve error checking
if exitFlag ~=1 && exitFlag ~= 2
%     warning('Possible problem solving system. Help!')
4;
end



end