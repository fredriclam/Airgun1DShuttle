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
%
% In final code:
% -- Check venting valve extreme cases (closed chamber at ambient pressure,
%    versus fully vented fixed ambient pressure)
% -- Orifice to middle chamber: check dimensions (area smaller than gap
%    characteristic size allows neglecting due to choked flow)