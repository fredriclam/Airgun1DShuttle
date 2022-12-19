function [agState, exception] = ...
    fullState(obj, q, t, bubble, shuttle, REVERT_MODEL, INCLUDE_ALL_PRIMITIVES)
% Computes the full state of the airgun from the state
% vector. Returns a named struct with subsystem variables.
%
% Use within RHS function for ODE and for querying state in
% post-processing.

exception = [];

if nargin <= 6
    % Disable to improve speed; enable to debug Euler domain states
    INCLUDE_ALL_PRIMITIVES = false;
end

% Tolerance for iterateToTol
iterativeSolveTol = 1e-6;
% Max number of iterations allowed for iterateToTol
iterateToTolMaxIterations = 1;
caseReiterationMax = 100;

if numel(q) > length(q)
    warning("fullState called for matrix of state variables. Function " + ...
            "is only intended for vector of state variables at a single t. ")
end

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

% Compute supersonic branch
function pRatio = pMinSonic_safe(gamma, ARatio)
    % Compute minimum pressure ratio p / pBubble for supersonic flow choked
    % at A_cs
    if gamma == 1.4
        % Map equation for pressure ratio to a polynomial equation, solve
        rr = roots([ARatio^-2,0,0,0,0,-(1+2/(gamma-1)),2/(gamma-1)]);
        supersonicRoot = rr(5);
        assert(isreal(supersonicRoot));
        assert(supersonicRoot >= 1 || ARatio < 1)
        % Map back to pressure ratio p / pBubble
        pRatio = supersonicRoot.^(7/2);
    else
        error("pMinSonic only configured for gamma = 1.4");
    end
end
pMinSonic = @(ARatio) pMinSonic_safe(obj.physConst.gamma, ARatio);

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
         warning('Entropy decreased from port to bubble.')
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
    [qPort, exitFlag] = ...
        obj.enforceScalarConstraint(essentialConstraint, qIn);
end

function qHat = processPortChokedCase()
    gamma = obj.physConst.gamma;
    % Map MPort = 1 -> A_port/A* -> MTarget
    [MTarget, exception]...
        = mapChokedPortAreaRatioToM(obj, APortExposed);

    numerator = u_R + 2/(gamma-1)  * c_R;
    denominator = MTarget * c_R + 2/(gamma-1) * c_R;
    cHat = numerator / denominator * c_R;
    uHat = MTarget * cHat;
    rhoHat = rho_R * (cHat / c_R)^(2/(gamma-1));
    THat = cHat^2 / gamma / obj.physConst.Q;
    eHat = (rhoHat * obj.physConst.c_v * THat + 0.5*rhoHat * uHat^2);
    qHat = [rhoHat;
            rhoHat * uHat;
            eHat];
end

function [qPort, exitFlag] = processChamberChokedCase(qIn)
    % True mach condition M_R = 1
    % This condition exists when there is lack of continuity between the
    % port-choked boundary case and the chamber-choked boundary case.
    % In the limit as dx -> 0 in the grid this case should disappear.
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
% Check the sonic condition
if pSonicFn(qIn) < pBubble
    caseKeyOut = 'subsonic';
    qOut = iterateToTol(...
        @processSubsonicCase, ...
        qIn, ...
        iterativeSolveTol);
else
    caseKeyOut = 'portChoked';
    qOut = processPortChokedCase();
    if qOut(3) < 0 || ...
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
% Check the supersonic branch minimum pressure

if obj.schm.p(qIn) / pBubble < ...
        pMinSonic(APortExposed / obj.physConst.crossSectionalArea)
    caseKeyOut = 'subsonic';
    qOut = iterateToTol(...
        @processSubsonicCase, ...
        qIn, ...
        iterativeSolveTol);
else
    % The following code block is an alternative strategy for enforcing sonic
    % boundary conditions (although it doesn't seem to produce physical
    % results)
    if false
        qOut = qIn;
        % Check grid value
        if M_R < 1
            caseKeyOut = 'subsonic';
            qOut = iterateToTol(...
                @processSubsonicCase, ...
                qIn, ...
                iterativeSolveTol);
        else
            caseKeyOut = 'chamberChoked';
        end
    end
    
    % Issue:
    % Does enforcing M == 1 eliminate the possibility of wave reflections
    % perturbing the boundary state from M = 1?
    if true
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
end

%% Compute boundary details
if REVERT_MODEL
    caseKey = 'noShuttle';
    % Fix outlet area to equal the cross-sectional area
    APortExposed = obj.physConst.crossSectionalArea;
    qPort = q_R;
    % Shutoff state:
    if t > obj.physConst.airgunCutoffTime
        uPort = qPort(2)/qPort(1);
        qPort(3) = qPort(3) - 0.5 * qPort(1) * uPort^2;
        qPort(2) = 0;
    end
else
    % Approximate the total port length as the full travel of
    % the shuttle: the % of the travel is thus the % of the
    % full port area that is exposed
%     APortExposed = max([0, obj.physConst.APortTotal * ...
%         (shuttle(1) - obj.physConst.portLead) / ...
%         (obj.physConst.operatingChamberLength - obj.physConst.portLead)]);

    % Appoximate exposed port area as
    %   xi - portLead
    % times the effective outer diameter
    obj.physConst.airgunPortLength = 2.375;
    APortExposed = (shuttle(1) - obj.physConst.portLead) ...
        / (0.0254*obj.physConst.airgunPortLength) ...
        * obj.physConst.APortTotal;
    % Clamp to [0, APortTotal]
    APortExposed = min([max([0, APortExposed]), obj.physConst.APortTotal]);
    F = 0.75;
    if isfield(obj.bubbleModel, 'flowreductionfactor')
        F = obj.bubbleModel.flowreductionfactor;
    end
    APortExposed = APortExposed * F;
    
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
    if APortExposed <= 1e-4 || t <= 0 % port is closed, or close enough
        % In case of t < 0, assume approximate equilibrium between gas in
        % and leakage.
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
                if reiterationCount > caseReiterationMax
                    error('Error finding case.')
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
                    error('Error finding case.')
                end
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

%% Post-process state
% Compute port sound speed
cPort = sqrt(obj.physConst.gamma * obj.physConst.Q * TPort);

% Bubble state recompute
% This calls bubbleRHS once the port state (in the form of a hat-state)
% is known. This is used to compute energy transfer terms using bubbleRHS.
[~, dQdt, workrate, dEin]= bubbleRHS(t, bubble, rhoPort, velocityPort, ...
    eTotalPort, pPort, ...
    obj.physConst.crossSectionalArea, obj.physConst, obj.bubbleModel);

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

%% Compute state at the port
% 

if ~REVERT_MODEL
    if strcmpi('portClosed', caseKey)
        p_outlet = pSonicFn(qPort);
    else
        Asonic_ = obj.physConst.crossSectionalArea / ...
            areaMachFunction(obj.physConst.gamma, MPort);
        if strcmpi('chamberChokedForced', caseKey)
            gamma = obj.physConst.gamma;
            M_outlet = 1;
%             M_outlet = fzero( @(M) ...
%                 ((gamma+1)/2)^(-(gamma+1)/2/(gamma-1)) * ...
%                 (1 + (gamma-1)/2 * M^2 )^ ...
%                 ((gamma+1)/2/(gamma-1)) ./ M - ...
%                 APortExposed / Asonic_, ...
%                 [1+1e-14,50]);
        else
            M_outlet = obj.machAreaFunction(APortExposed / Asonic_);
        end
        pStagnation = pPort / pressureMachFunction(obj.physConst.gamma, MPort);
        p_outlet = pStagnation * pressureMachFunction(obj.physConst.gamma, M_outlet);
    end
else
    p_outlet = pPort;
end

%% Shuttle detailed state
if length(shuttle) >= 2
[~, ~, ~, ~, subsystemState] = shuttleEvolve(shuttle, ...
    pPort, obj.physConst, obj.chambers, rhoPort * cPort);
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
    'pSonicPort', pSonicPort, ...
    'p_outlet', p_outlet ...
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
    'rho', rhoBubble, ...
    'dQdt', dQdt, ...
    'pdV', workrate, ...
    'dEin', dEin ...
);
bubbleStates.rho = bubbleStates.p / ...
    (obj.physConst.Q * bubbleStates.T);
if length(shuttle) >= 2
    shuttleStates = subsystemState;
else
    shuttleStates = struct();
end

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

