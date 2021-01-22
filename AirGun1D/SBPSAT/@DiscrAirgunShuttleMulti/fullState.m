function [agState, exception] = ...
    fullState(obj, q, t, bubble, shuttle, REVERT_MODEL)
% Computes the full state of the airgun from the state
% vector. Returns a named struct with subsystem variables.
%
% Use within RHS function for ODE and for querying state in
% post-processing.

exception = [];
INCLUDE_ALL_PRIMITIVES = true;
USE_LEGACY_ENFORCEMENT = false;

% Compute primitive variables at right of PDE domain
q_R = obj.schm.e_R'*q;
p_R = obj.schm.p(q_R);
u_R = q_R(2)/q_R(1);
rho_R = q_R(1);
e_tot_R = q_R(3);
T_R = (e_tot_R - 0.5 * rho_R * u_R^2) / rho_R / obj.physConst.c_v;
c_R = sqrt(obj.physConst.gamma * obj.physConst.Q * T_R);
M_R = u_R / c_R;

% Mach factor
machFactor = @(M) (1 + (obj.physConst.gamma - 1)/2 * ...
    M^2);
% Compute stagnation pressure
p0_R = p_R * ...
    (machFactor(M_R))^(obj.physConst.gamma/(obj.physConst.gamma-1));
% Compute sonic pressure
pSonic_R = p0_R * ...
    (machFactor(1))^(-obj.physConst.gamma/(obj.physConst.gamma-1));

% Capture flow state in PDE domain
flowState = obj.schm.flowStateR(q);

% Compute bubble density
pBubble = bubblePressure(bubble, obj.physConst);
TBubble = bubble(4) / obj.physConst.c_v / bubble(3);
rhoBubble = pBubble / (obj.physConst.Q * TBubble);

% Create empty note string
noteString = "";

% TEMP TODO: replace
qPort = q_R;
% TEMP TODO: default
% Compute port characteristic variable values
wPort = obj.schm.T(qPort) \ qPort;

% Compute geometry
if REVERT_MODEL
    caseKey = 'noShuttle';
    % Fix outlet area to equal the cross-sectional area
    APortExposed = obj.physConst.A;
    if t <= obj.physConst.AirgunCutoffTime
        massFlowPort = rho_R * u_R * APortExposed;
        velocityPort = u_R;
    else
        massFlowPort = 0;
        velocityPort = 0;
    end
    rhoPort = rho_R;
    pPort = p_R;
    TPort = T_R;
else
    % Approximate the total port length as the full travel of
    % the shuttle: the % of the travel is thus the % of the
    % full port area that is exposed
    APortExposed = max([0, obj.physConst.APortTotal * ...
        (shuttle(1) - obj.physConst.portLead) / ...
        (obj.physConst.operatingChamberLength - obj.physConst.portLead)]);
    
    % Case dependent boundary values
    if APortExposed == 0 % Case 0A: port is closed
        caseKey = 'portClosed'; % Port is closed
        massFlowPort = 0;
        velocityPort = 0;
        pPort = p_R;
        TPort = T_R;
        rhoPort = pPort / obj.physConst.Q / TPort;
    elseif obj.schm.flowStateR(q) == scheme.Euler1d.SUPERSONIC_OUTFLOW
                                     % Case 0B: natural (super)sonic outflow
        caseKey = 'chamberChokedNatural';
        
        velocityPort = u_R;
        pPort = p_R;
        TPort = T_R;
        rhoPort = pPort / obj.physConst.Q / TPort;
        massFlowPort = rhoPort * velocityPort * obj.physConst.crossSectionalArea;
    elseif pSonic_R < pBubble % Case 1: insufficient pressure to choke at port
        caseKey = 'subsonic';
        % TODO: Replace with simultaneous characteristic enforcement
        rhoPort = rho_R;  % Not quite right
        TPort = T_R;      % Not quite right
        % Compute upstream pressure
        pUpstream = @(pDownstream, rhoDownstream, rhoUpstream) ...
          pDownstream * (rhoUpstream / rhoDownstream)^obj.physConst.gamma;
        pPort = pUpstream(pBubble, rhoBubble, rhoPort);
        % Compute upstream mach from pressure
        try
            MPort = machPressureFunction(obj.physConst.gamma, pPort/p0_R);
        catch
            caseKey = 'subsonicForcingSonic';
            % DOCUMENT: Use of simple enforcement is not quite right and
            % leads to these jump cases
            MPort = 1;
        end
        velocityPort = MPort * c_R;
        massFlowPort = rho_R * velocityPort * ...
            obj.physConst.crossSectionalArea;
    elseif APortExposed <= obj.physConst.crossSectionalArea
        caseKey = 'portChoked';
        [MPort, exception]...
            = mapChokedPortAreaRatioToM(obj, APortExposed);
        
        if USE_LEGACY_ENFORCEMENT
            velocityPort = MPort * c_R;
            massFlowPort = velocityPort * rho_R * obj.physConst.crossSectionalArea;
            rhoPort = rho_R;
            pPort = p_R;
            TPort = T_R;
        else
            % Create namespace
            enforcement = struct();
            % Projection matrix to outgoing characterisics
            enforcement.P = [1 0 0
                             0 1 0];
            enforcement.mapq2M = @(q) sign(q(2)) * sqrt(...
                (1/((obj.physConst.gamma - 1)*obj.physConst.gamma))...
                * (q(2)^2) / (q(1)*q(3) - 0.5*q(2)^2) ...
            );
            enforcement.mapq2characteristics = @(q) obj.schm.T(q) \ q;
            
            essentialConstraint = @(q) ...
                enforcement.mapq2M(q) - MPort;
            outgoingCharBCConstraints = @(q) ...
                enforcement.P ...
                * (enforcement.mapq2characteristics(q) ...
                   - enforcement.mapq2characteristics(q_R));
            % Add scaling to equation system to give weight to essential
            % constraint
            scalingMatrix = diag([max(1,min(1e5,1/M_R)), 1, 1]);
            scaledConstraintSystem = @(q) scalingMatrix * [
                essentialConstraint(q);
                outgoingCharBCConstraints(q);
            ];
        
            % Solve constraint system
            [qPort, constraintVals, exitFlag] = fsolve(...
                scaledConstraintSystem, q_R, ...
                optimoptions('fsolve','display','none'));
            if exitFlag ~=1
                warning('Possible problem solving system. Help!')
            end
            
            % Compute primitives
            pPort = obj.schm.p(qPort);
            velocityPort = qPort(2)/qPort(1);
            rhoPort = qPort(1);
            eTotalPort = qPort(3);
            TPort = (eTotalPort - 0.5 * rhoPort * velocityPort^2) / rhoPort / ...
                obj.physConst.c_v;
            cPort = sqrt(obj.physConst.gamma * obj.physConst.Q * TPort);
            MPort = velocityPort / cPort;
            massFlowPort = rhoPort*velocityPort*obj.physConst.crossSectionalArea;
            wPort = enforcement.mapq2characteristics(qPort);
%             disp(wPort(3));
        end
    elseif APortExposed > obj.physConst.crossSectionalArea
        caseKey = 'chamberChokedForced';
        % TODO: Replace with properly enforced Mach condition
        MPort = 1;
        velocityPort = MPort * c_R;
        massFlowPort = velocityPort * rho_R * obj.physConst.crossSectionalArea;
        rhoPort = rho_R;
        pPort = p_R;
        TPort = T_R;
    end
    
%     if  flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW ...
%             && ~strcmpi ('chamberChoked', caseKey)
%         noteString = noteString + ...
%             "Chamber should not be choked but it is";
%     end
    
    % Compute port total energy per unit volume
    eTotalPort = rhoPort * obj.physConst.c_v * TPort ...
        + 0.5 * rhoPort * velocityPort^2;
end

% Compute port sound speed
cPort = sqrt(obj.physConst.gamma * obj.physConst.Q * TPort);

if INCLUDE_ALL_PRIMITIVES
    % Compute primitives everywhere in Euler domain
    pEulerDomain = NaN(size(q(1:3:end)));
    for i = 1:length(pEulerDomain)
        pEulerDomain(i) = obj.schm.p(q(1+3*(i-1):3+3*(i-1)));
    end
    rhoEulerDomain = q(1:3:end);
    uEulerDomain = q(2:3:end) ./ q(1:3:end);
    eEulerDomain = q(3:3:end);
    TEulerDomain = (eEulerDomain ./ rhoEulerDomain ...
        - 0.5 * uEulerDomain.^2) / obj.physConst.c_v;
    cEulerDomain = sqrt(obj.physConst.gamma * obj.physConst.Q * TEulerDomain);
    MEulerDomain = uEulerDomain ./ cEulerDomain;
else
    pEulerDomain = [];
    rhoEulerDomain = [];
    uEulerDomain = [];
    eEulerDomain = [];
    TEulerDomain = [];
    cEulerDomain = [];
    MEulerDomain = [];
end

%% Construct data hierarchy
portStates = struct( ...
    'APortExposed', APortExposed, ...
    'massFlowPort', massFlowPort, ...
    'velocityPort', velocityPort, ...
    'rhoPort', rhoPort, ...
    'pPort', pPort, ...
    'TPort', TPort, ...
    'cPort', cPort, ...
    'MPort', velocityPort/cPort, ...
    'eTotalPort', eTotalPort, ...
    'caseKey', caseKey, ...
    'ARatio', APortExposed / obj.physConst.crossSectionalArea, ...
    'qPort', qPort, ...
    'wPort', wPort ...
);
eulerDomainStates = struct(...
    'rho', rhoEulerDomain, ...
    'u', uEulerDomain, ...
    'e', eEulerDomain, ...
    'p', pEulerDomain, ...
    'T', TEulerDomain, ...
    'c', cEulerDomain, ...
    'M', MEulerDomain, ...
    'q_R', q_R, ...
    'p_R', p_R, ...
    'u_R', u_R, ...
    'rho_R', rho_R , ...
    'e_tot_R', e_tot_R, ...
    'T_R', T_R, ...
    'c_R', c_R, ...
    'M_R', M_R, ...
    'p0_R', p0_R, ...         % Stagnation pressure
    'pSonic_R', pSonic_R, ... % Sonic pressure
    'flowState', flowState, ...
    'matrix_T_at_x_R', obj.schm.T(q_R), ...
    'characteristics_at_x_R', obj.schm.T(q_R) \ q_R ...
);
bubbleStates = struct(...
    'R', bubble(1), ...
    'RDot', bubble(2), ...
    'm', bubble(3), ...
    'E', bubble(4) , ...
    'p',  pBubble, ...
    'T', TBubble, ...
    'rho', rhoBubble ...
);
bubbleStates.rho = bubbleStates.p / ...
    (obj.physConst.Q * bubbleStates.T);
shuttleStates = struct(...
    'pos', shuttle(1), ...
    'vel', shuttle(2) ...
);

agState = struct(...
    't', t, ...
    'portStates', portStates, ...
    'eulerDomainStates', eulerDomainStates, ...
    'bubbleStates', bubbleStates, ...
    'shuttleStates', shuttleStates, ...
    'noteString', noteString ...
);
disp(caseKey);
    
%     figure(103);
%     bar(eulerDomainStates.characteristics_at_x_R)
%     set(gca, 'XTickLabel', {'w_{u}','w_{u+c}','w_{u-c}'})
%     drawnow
end