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