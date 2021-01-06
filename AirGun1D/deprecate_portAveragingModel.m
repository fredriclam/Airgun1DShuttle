% A possible model for modeling the port area, where an averaged state
% over the isentropic contraction is computed to estimate the pressure
% acting on the shuttle.

% If bubble pressure > port pressure, the
% flow is not consistent; switch to subsonic.
% TODO: check code for subsonic case
if p_R < pBubble
    warning(['p_R < pBubble. Should do subsonic flow; t = ' num2str(t)]);
    
    % Inverse problem for finding subsonic mach number at port
    obj_func = @(M_port) APortExposed / ...
        physConst.cross_sectional_area- ...
        area_ratio(M_port) / area_ratio(M_R);
    % Special case: M is zero
    if ~isfinite(area_ratio(M_R))
        M_port = 0;
    else
        % Solve for subsonic mach root
        M_port = fzero(obj_func, [1e-10,1 - 1e-10]);
        assert(M_port >= 0 && M_port <= 1);
    end
    % Compute thermodynamic properties using forward
    % relations
    pPort = pressure_ratio(M_port) / ...
        pressure_ratio(M_R) * p_R;
    % Compute pressure upstream through isentropic
    % flow relation
    p_a = pressure_ratio(M_R) / ...
        pressure_ratio(M_port) * pBubble;
    % Compute temperature downstream through isentropic
    % flow relation
    TPort = temperature_ratio(M_port) / ...
        temperature_ratio(M_R) * T_R;
    rhoPort = pPort / physConst.Q / TPort;
    % Averaged density over the control volume
    % (substitute)
    rhoPortAveraged = mean([rho_R, rhoPort]);
    % Suppose the mass flow rate is determined by inflow to
    % control volume minus volumetric dilation; use
    % averaged density between two endpoints of isentropic
    % flow as approximation to the bulk-averaged density
    massFlowPort = physConst.cross_sectional_area * ...
        (- rhoPort * velShuttle ...
        ...(- rhoPortAveraged * velShuttle ...
        + rho_R * u_R);
    % Explicit subsonic
    isSonicFlags(2) = false;
else
    % Flag port as sonic
    isSonicFlags(2) = true;
end