%% Notes
%
% Updated shuttleEvolve.m: include human-readable data struct for proper
%   state monitoring
%
% To do:
% -- Check and writeup formula for flow within operating chamber
%      shuttleEvolve.m >> line:
%           flowL2R = ...
% -- Monitor exports from ODE solver (or dummy state) instead of using RHS
%      raw output
% -- Pass up the shuttle stuff
% -- Clean up function to reset the bubble frozen boolean
% -- Check the flow work in both ways for the op chamber
%
% -- Animation to show the BC switching
% --
% -- Note: assume that the port ramps up from 0% to 100% over the entire
%      operating chamber length. So the operating chamber length is a bit
%      shorter than it is. Let's check this.
%
% In final code:
% -- Check venting valve extreme cases (closed chamber at ambient pressure,
%    versus fully vented fixed ambient pressure)
% -- Orifice to middle chamber: check dimensions (area smaller than gap
%    characteristic size allows neglecting due to choked flow)
% -- Check sensitivty to the threshold opening distance. This seems rather
%    important to the final equilibrium position.

% Issues:
% -- Enforcement with characteristics (code is built for only incoming
% characteristic enforcement) makes the delta w << w (because delta q << q)
% so there may be a scaling issue
% -- Interpretation of the system-solve as a "locally compatible" boundary
% state allows us to place any boundary type consistent with the
% system-solve (I think). Let's see what we think about this.
% -- Some contributors to negative u:
%    -- Initial wigglies for double-wall condition (left and right closed)
%    -- Tiny drop in pressure, w required for the initial acceleration of
%    gas. Need to enforce in u, and look into this further.

% Notes:
% -- Sensitivity of w solve:
% Scaling matrix of [M - MPort; w_u - w_uPort; w_{u+c} - w_{u+c}Port] varies
% in the following (same tolerances)
%
%   Scaling values of [1e5, 1, 1] lead to q of
%   1.0e+07 *
%   0.000001507752336
%   0.000016481027984
%  1.724523832767844
% 
%  qStar1 = fsolve(...
%    @(q) diag([1 1 1]) * [
% essentialConstraint(q);
% outgoingCharBCConstraints(q);
% ], q_R, optimoptions('fsolve','FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13));
%  
% 
% qStar1 =
% 
%    1.0e+07 *
% 
%    0.000001507679872
%    0.000016480573451
%    1.724511594529822
% 
%  qStar2 = fsolve(...
% @(q) diag([1 1e5 1e5]) * [
% essentialConstraint(q);
% outgoingCharBCConstraints(q);
% ], q_R, optimoptions('fsolve','FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13));
% 
% 
% Equation solved, solver stalled.
% 1.0e+07 *
% 
%    0.000003856089409
%    0.000022398504588
%    1.245436385885649
