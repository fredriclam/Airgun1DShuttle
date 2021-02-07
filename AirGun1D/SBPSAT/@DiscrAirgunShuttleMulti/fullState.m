function [agState, exception] = ...
    fullState(obj, q, t, bubble, shuttle, REVERT_MODEL)
% Computes the full state of the airgun from the state
% vector. Returns a named struct with subsystem variables.
%
% Use within RHS function for ODE and for querying state in
% post-processing.

% TODO: track #evals of each case

exception = [];
% Disable to improve speed; enable to debug Euler domain states
INCLUDE_ALL_PRIMITIVES = false;

% Tolerance for iterateToTol
iterativeSolveTol = 1e-8;
% Max number of iterations allowed for iterateToTol
iterateToTolMaxIterations = 100;

%% Compute primitive variables at right of PDE domain
q_R = obj.schm.e_R'*q;
p_R = obj.schm.p(q_R);
u_R = q_R(2)/q_R(1);
rho_R = q_R(1);
e_tot_R = q_R(3);
T_R = (e_tot_R - 0.5 * rho_R * u_R^2) / rho_R / obj.physConst.c_v;
c_R = obj.schm.c(q_R);
M_R = u_R / c_R;

%% Compute isentropic flow quantities
% Mach factor
machFactor = @(M) (1 + (obj.physConst.gamma - 1)/2 * ...
    M^2);
% Mach number function (signed)
machFn = @(q) (q(2)/q(1)) / obj.schm.c(q);
% Compute stagnation pressure
pStagnationFn = @(q) obj.schm.p(q) * ...
    machFactor(machFn(q))^(obj.physConst.gamma/(obj.physConst.gamma-1));
p0_R = p_R * ...
    (machFactor(M_R))^(obj.physConst.gamma/(obj.physConst.gamma-1));
% Compute sonic pressure
pSonicFn = @(q) pStagnationFn (q) * ...
    (machFactor(1))^(-obj.physConst.gamma/(obj.physConst.gamma-1));
pSonic_R = pSonicFn(q_R);
% pSonic_R = p0_R * ...
%     (machFactor(1))^(-obj.physConst.gamma/(obj.physConst.gamma-1));

%% Capture flow state in PDE domain
flowState = obj.schm.flowStateR(q);

%% Compute bubble state
pBubble = bubblePressure(bubble, obj.physConst);
TBubble = bubble(4) / obj.physConst.c_v / bubble(3);
rhoBubble = pBubble / (obj.physConst.Q * TBubble);

%% Define map functions
% Define map to Mach number (signed)
mapq2M = @(q) sign(q(2)) * sqrt(...
    (1/((obj.physConst.gamma - 1)*obj.physConst.gamma))...
    * (q(2)^2) / (q(1)*q(3) - 0.5*q(2)^2) ...
);
% Define map to approximate characteristics
mapq2characteristics = @(q) obj.schm.T(q) \ q;

%% Create empty note string
noteString = "";

% % Default value
% qPort = q_R;

%% Define function in scope of fullState
% Functions that handle boundary cases, capturing the local workspace
function [qPort, exitFlag] = processSubsonicCase(qIn)
    % Note that entropy is lower in the outlet flow than in the bubble;
    % mixing, shocks, and turbulent dissipation causes a subsequent
    % entropy increase that is not explicitly modeled.
    % Thus s_R <= sPort <= sBubble

    % Check entropy (exponential-entropy)
    entropyFn = @(p, rho) p / rho^obj.physConst.gamma;
    if entropyFn(p_R, rho_R) >= entropyFn(pBubble, rhoBubble)
        error('Entropy decreased from port to bubble.')
    end

    % Build the function for "Mach number consistent with
    % downstream pressure continuity"
    % Map from q to mach at the port
    machExitFn = @(q) machPressureFunction(obj.physConst.gamma, ...
          pBubble/pStagnationFn(q));
    % Composed map: q -> mach at the port -> A_port/A* -> A_cs/A* ->
    %   M_upstream
    MPortFn = @(q) obj.machAreaFunction(...
        obj.physConst.crossSectionalArea / ...
        APortExposed * areaMachFunction(obj.physConst.gamma, ...
        machExitFn(q)));

    % Constraint: mach computed from q is equal to the mach computed
    % from steady-state flow terminating in pressure continuity at the
    % port
    essentialConstraint = @(q) machFn(q) - MPortFn(q);

    % Solve for viable q
    % TODO: prove uniqueness?
    [qPort, exitFlag] = ...
        obj.enforceScalarConstraint(essentialConstraint, qIn);
end

function [qPort, exitFlag] = processPortChokedCase(qIn)
    [MPort, exception]...
        = mapChokedPortAreaRatioToM(obj, APortExposed);

    % Mach boundary condition at exit of PDE domain
    essentialConstraint = @(q) ...
            mapq2M(q) - MPort;

    [qPort, exitFlag] = ...
        obj.enforceScalarConstraint(essentialConstraint, qIn);
end

function [qPort, exitFlag] = processChamberChokedCase(qIn)
    % True mach condition M_R = 1
    % This condition exists if the Euler domain feels it should be
    % subsonic, but such a setup would result in a subsonic expansion
    % into the bubble, where the sonic pressure is too high and the
    % flow would have choked--precisely at the exit of the firing
    % chamber. This makes sure M_R < 1 is invalid, and corrects the
    % PDE. This boundary may or may not be represented in the final ode
    % system solution, but has been encountered during ode solve.
    essentialConstraint = @(q) ...
            mapq2M(q) - 1;

    [qPort, exitFlag] = ...
        obj.enforceScalarConstraint(essentialConstraint, qIn);
end

function [qNew, exitFlag, iterateCount] = iterateToTol(callback, q, tol)
    % Iterate a contraction mapping to tolerance in the norm
    iterateCount = 1;
    [qNew, exitFlag] = callback(q);
    if exitFlag ~= 1
        return
    end
    
    while norm(qNew-q) > tol && iterateCount < iterateToTolMaxIterations
        q = qNew;
        [qNew, exitFlag] = callback(q);
        if exitFlag ~= 1
            return
        end
        iterateCount = iterateCount + 1;
    end
    
%     if iterateCount >= iterateToTolMaxIterations
%         warning('Iteration count reached in BC inner loop')
%     end
end

% Select subsonic or sonic case for contraction flow (cross section area >
% port area).
function [qOut, caseKeyOut] = tryCasesContraction(qIn)
if pSonicFn(qIn) < pBubble
    caseKeyOut = 'subsonic';
    qOut = iterateToTol(...
        @processSubsonicCase, ...
        qIn, ...
        iterativeSolveTol);
else
    caseKeyOut = 'portChoked';
    [qOut, exitFlag, ~] = iterateToTol(...
        @processPortChokedCase, ...
        qIn, ...
        iterativeSolveTol);
    if exitFlag ~= 1 || exitFlag == -9 || ...
            qOut(3) < 0 || ...
            ~(obj.schm.c(qOut) > 0) || ...
            ~(obj.schm.p(qOut) > 0)
        % Relax the numerics (required for predictor steps
        % sometimes)
        % Issue:
        % M_R >> 1 for a step during a portClosed-portChoked only
        % sequence. See commit [wip 4b534a7].
        % Fallback case when M_R is too far away from
        % subsonic, and we fail to find a q such that M(q) = 1
        % while preserving the outgoing characteristics.
        caseKeyOut = 'relaxation';
        qOut = qIn;
    end
end
end

% Select subsonic or sonic case for expansion flow (cross section area >
% port area).
function [qOut, caseKeyOut] = tryCasesExpansion(qIn)
if pSonicFn(qIn) < pBubble
    caseKeyOut = 'subsonic';
    qOut = iterateToTol(...
        @processSubsonicCase, ...
        qIn, ...
        iterativeSolveTol);
else
    caseKeyOut = 'chamberChokedForced';
    % Mach-1 boundary condition at exit of PDE domain
    [qOut, exitFlag, ~] = iterateToTol(...
        @processChamberChokedCase, ...
        qIn, ...
        iterativeSolveTol);
    if exitFlag ~= 1 || ...
            qOut(3) < 0 || ...
            ~(obj.schm.c(qOut) > 0) || ...
            ~(obj.schm.p(qOut) > 0)
        % Relax the numerics (required for predictor steps
        % sometimes)
        % Issue:
        % M_R >> 1 for a step during a portClosed-portChoked only
        % sequence. See commit [wip 4b534a7].
        % Fallback case when M_R is too far away from
        % subsonic, and we fail to find a q such that M(q) = 1
        % while preserving the outgoing characteristics.
        caseKeyOut = 'relaxation';
        qOut = qIn;
    end
end
end

%% Compute boundary details
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
    
    % Compute initial q from q_R such that M_R as <= 1
    if M_R > 1
        % Iteratively solve for q
        [q_RMod, exitFlag] = iterateToTol(...
            @processChamberChokedCase, ...
            q_R, ...
            iterativeSolveTol);
        % Replace q_R with modified q if such a solution is found
        if exitFlag == 1
            q_R = q_RMod;
        end
    end
    
    % Compute new q_R depending on boundary case
    if APortExposed == 0 % port is closed
        % Kinematically-constrained case; no boundary values necessary
        % (wall closure used)
        caseKey = 'portClosed'; % Port is closed
        % Compute q-hat such that velocity is removed
        % Note: this value is not used for boundary enforcement (see
        % closure function instead)
        qPort = [rho_R;
                 0;
                 e_tot_R - 0.5 * rho_R * u_R^2];
    else
        if APortExposed <= obj.physConst.crossSectionalArea % (contraction)
            % Select sonic or subsonic treatment based on q_R
            [qPort, caseKey] = tryCasesContraction(q_R);
            
            caseKeyOld = '';
            % Outer iteration loop: until the boundary case stabilizes
            reiterationCount = 1;
            while ~strcmpi(caseKeyOld, caseKey)
                % Swap data
                caseKeyOld = caseKey;
                % Select sonic or subsonic treatment based on qPort
                [qPort, caseKey] = tryCasesContraction(qPort);
                
                reiterationCount = reiterationCount + 1;                
                if reiterationCount > 100
                    error('Cant decide which case this is. Help!')
                end
            end
        else % APortExposed > obj.physConst.crossSectionalArea (expansion)
            % Select sonic or subsonic treatment based on q_R
            [qPort, caseKey] = tryCasesExpansion(q_R);
            
            caseKeyOld = '';
            % Outer iteration loop: until the boundary case stabilizes
            reiterationCount = 1;
            while strcmpi(caseKeyOld, caseKey)
                % Swap data
                caseKeyOld = caseKey;
                % Select sonic or subsonic treatment based on qPort
                [qPort, caseKey] = tryCasesExpansion(qPort);
                
                reiterationCount = reiterationCount + 1;                
                if reiterationCount > 100
                    error('Cant decide which case this is. Help!')
                end
            end
        end
    end
    
    % Compute primitives
    rhoPort = qPort(1);
    velocityPort = qPort(2) / qPort(1);
    eTotalPort = qPort(3);
    TPort = (eTotalPort - 0.5 * rhoPort * velocityPort^2) / rhoPort / ...
        obj.physConst.c_v;
    pPort = obj.schm.p(qPort);
    cPort = obj.schm.c(qPort);
    MPort = velocityPort / cPort;
    massFlowPort = rhoPort*velocityPort*obj.physConst.crossSectionalArea;
    wPort = mapq2characteristics(qPort);
    pSonicPort = pSonicFn(qPort);
end

%% Post-process state
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
    'wPort', wPort, ...
    'pSonicPort', pSonicPort ...
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
% disp(caseKey);
    
%     figure(103);
%     bar(eulerDomainStates.characteristics_at_x_R)
%     set(gca, 'XTickLabel', {'w_{u}','w_{u+c}','w_{u-c}'})
%     drawnow
end

