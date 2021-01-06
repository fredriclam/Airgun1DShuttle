% Compute port dynamics based on mass conservation on control
% volume (0-dimensional model with average rho or T dependent on
% upstream or downstream parameters).
function massConservationBasedPhysics()
error('Don''t use me. Not fully implemented')
if false % Mass conservation choking
    % Compute subsonic mass flow rate
    pPort = pBubble; % Match downstream
    rhoPort = pPort / physConst.Q / TPort;
    massFlowPort = physConst.cross_sectional_area * ...
        (- rhoPort * velShuttle ...
        + rho_R * u_R);
    
    % IDEA: replace u_R with near-boundary average and
    % smoothing
    uPortGuess = massFlowPort / rhoPort / APortExposed;
    c_port = sqrt(physConst.gamma * physConst.Q * ...
        TPort);
    %                 % Messwithit
    % %                 rhoPort = pBubble / physConst.Q / TBubble;
    %                 uPortGuess = massFlowPort / rhoPort / APortExposed;
    %                 c_port = sqrt(physConst.gamma * physConst.Q * ...
    %                     TBubble);
    % Enforce port sonic -> subsonic transition once only
    if t < 0 % 1e-6 % "subsonic" only
        isSonicFlags(2) = false;
    else
        if miscStates(1) == 0 % Accelerating flow
            if uPortGuess >= c_port
                % True sonic regime
                miscStates(1) = 1;
            end
            isSonicFlags(2) = false;
        elseif miscStates(1) == 1 % Stabilization-needed state
            %                     % Can transition to subsonic only
            if uPortGuess < 0.99*c_port && uPortGuess > 0.0*c_port
                % True sonic regime
                miscStates(1) = 2;  %miscStates(1) = 1;
            else
                isSonicFlags(2) = true;
            end
            %                         isSonicFlags(2) = true;
        end
        
        if miscStates(1) == 2 % Subsonic flow
            if uPortGuess >= 1.01*c_port
                miscStates(1) = 1;
                isSonicFlags(2) = true;
            else
                isSonicFlags(2) = false;
            end
        end
    end
    
    % Compute rho * velocity in airgun
    rhovel_a = massFlowPort / ...
        physConst.cross_sectional_area + ...
        rhoPort * velShuttle;
end
end