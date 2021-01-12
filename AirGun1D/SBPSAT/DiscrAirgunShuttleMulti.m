classdef DiscrAirgunShuttleMulti < DiscrAirgun
    properties
        shuttle0            % Initial shuttle state [pos; vel]
        plug0               % Initial plug control volume state [mass; en]
        machAreaFunction    % Precomputed M(A/A*) function interpolant
        opChamber
        bubbleFrozen        % Freeze bubble state until first open to prevent negative bubble pressure
    end
    
    methods
        function obj = DiscrAirgunShuttleMulti(nx,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth, airgunCrossSectionalArea, REVERT_MODEL, ...
                airgunFiringChamberProfile, ...
                airgunOperatingChamberProfile, bubbleInitialVolume, ...
                shuttleBdryPenaltyStrength)
            if nargin == 8 && ~REVERT_MODEL
                error('Incorrect # of arguments for unreverted model.')
            end
            
            DEBUG = false;
            
            % Call parent constructor
            obj = obj@DiscrAirgun(nx,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            obj.opChamber = OperatingChamber();
            if ~REVERT_MODEL
                % Override the configuration with updated, backward-compatible
                % setting
                [physConst, t0, icAirgun, icBubble] = ...
                    configAirgun('GeneralAirgun', ...
                        airgunPressure, ...
                        airgunLength, ...
                        airgunPortArea, ...
                        airgunDepth, ...
                        airgunCrossSectionalArea, ...
                        airgunFiringChamberProfile, ...
                        airgunOperatingChamberProfile, ...
                        bubbleInitialVolume, ...
                        shuttleBdryPenaltyStrength);
                obj.physConst = physConst;
            end
            
            % Alias commonly used objects
            physConst = obj.physConst;
            schm = obj.schm;
            gamma_ = physConst.gamma;
            Q_ = physConst.Q;
            % Replace this object's description
            obj.description = ...
                'Airgun augmented with port-region CV and shuttle';
            
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);     % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0; % Temperature
            p0 = physConst.p0a;   % Initial pressure
            
            % Define initial shuttle state: position; velocity
            % NOTE: position non-zero or else singular--for CV treatment
            % TODO: check sensitivity
            % TODO: incorporate initial volume encompassed by shuttle but
            % not relevant to the port area calculation
            if REVERT_MODEL
                V_front_max = 0;
                r = 1;
                mass_operatChamberAir = 0;
                x0_rear = 0;
            else
                V_front_max = physConst.shuttle_area_right * ...
                              physConst.operatingChamberLength;
                x0_rear = 1e-3;
                mass_operatChamberAir = rho0*V_front_max;
                r = x0_rear / physConst.operatingChamberLength;
                disp(['Mass in rear partition: ' ...
                      num2str(r * mass_operatChamberAir)]);
            end
            
            obj.shuttle0 = [x0_rear; % [m] -- must give rear chamber some room
                            0;
                            r * mass_operatChamberAir;
                            physConst.c_v*T0*(r * mass_operatChamberAir);
                            (1-r)*mass_operatChamberAir;
                            physConst.c_v*T0*(1-r)*mass_operatChamberAir;
                            ];   % [m/s]

            if ~REVERT_MODEL
                obj.plug0 = [rho0*physConst.plugVolumeFn(obj.shuttle0(1)); % [kg]
                    rho0*physConst.plugVolumeFn(obj.shuttle0(1))* ...
                    physConst.c_v*T0];                      % [J]
            else
                obj.plug0 = [0; 0];
            end
            
            obj.bubbleFrozen = true;
                       
                     
            %% Create boundary condition operators
            closure_l = schm.boundary_condition('l', 'wall');
            % Outflow with pressure (for unchoked-everywhere flow)
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');
            % Outflow with velocity (for choked flow)
            closure_r_out_sub_vel = ...
                schm.boundary_condition('r', 'outflow_vel');
            % Outflow with rho*velocity (alternative to the latter; unused)
            closure_r_out_sub_rhovel = ...
                schm.boundary_condition('r', 'outflow_rhovel');
            % Wall on the right
            closure_r_closed = schm.boundary_condition('r', 'wall');
            % Precompute mach area function M(A/A*) and, for subsonic,
            % the mach pressure function M(p/p0)
            if ~REVERT_MODEL
                obj.machAreaFunction = precomputeMachAreaFunction(gamma_);
            end
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle, dPlug, p_RTarget, u_RTarget, monitor] = ...
                    RHS(q, t, bubble, shuttle, plug, p_RTarget, u_RTarget)
                %% Set flag for bypassing plug modeling
                BYPASS_PLUG =  true;
                USE_SHIFT = false;
                
                %% Precompute
                % Compute primitive variables at right of PDE domain
                if REVERT_MODEL
                    q_R = schm.e_R'*q;
                else
%                     q_R = q(end-5:end-3); % Go one node in for weak enforcement effects
%                     if sign(q(end-1)) * sign(q_R(2)) < 0
%                         q_R(2) = 0;
%                     end
%                     q_R = schm.e_R'*q;
                    q_R = q(end-5:end-3); % Go one node in for weak enforcement effects
                end

                % Replace p_R with target pressure
                if REVERT_MODEL
                    p_R = schm.p(q_R);
                    u_R =  q_R(2)/q_R(1);
                else
                    % p_RTarget empty means no pressure enforcement
                    if ~isempty(p_RTarget)
                        p_R = p_RTarget;
                    else
                        p_R = schm.p(q_R);
                        if p_R < 0
                            p_R = schm.p(q(end-8:end-6));
                        end
                    end
                    
                    if ~isempty(u_RTarget)
%                         u_R = u_RTarget; % TODO: USE
                        u_R = q_R(2)/q_R(1);
                    else
                        u_R = q_R(2)/q_R(1);
                    end
                    
                                            
                end
                % Default to empty target pressure, velocity constraint
                p_RTarget = [];
                u_RTarget = [];
                
                rho_R = q_R(1);
                en_R = q_R(3);
                T_R = p_R / rho_R / Q_;
                c_R = sqrt(physConst.gamma * Q_ * T_R);
                M_R = u_R / c_R;
                
%                 % "Fixing" the supersonic weak enforcement case
%                 if M_R > 1
%                     warning('Supersonic capped to sonic.')
%                     M_R = 1;
%                     u_R = M_R * c_R;
%                 end
                
                % Compute bubble variables
                pBubble = bubblePressure(bubble, physConst);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                rhoBubble = pBubble / Q_ / TBubble;
                % Extract shuttle variables
                posShuttle = shuttle(1);
                velShuttle = shuttle(2);
                % Compute geometry
                if REVERT_MODEL
                    % Fixed outlet area
%                     APortExposed = physConst.APortTotal;
                    % Fix outlet area to equal the cross-sectional area
                    APortExposed = physConst.A;
                else
                    % Approximate the total port length as the full travel of
                    % the shuttle: the % of the travel is thus the % of the
                    % full port area that is exposed
                    APortExposed = max([0, physConst.APortTotal * ...
                        (posShuttle - physConst.portLead) / ...
                        (physConst.operatingChamberLength - physConst.portLead)]);
%                     APortExposed = 0;

                    if USE_SHIFT
                        shiftCounter = 1;
                        while p_R < 0 || rho_R < 0 || T_R < 0 || M_R < -1e-2
                            warning('Used next point in boundary.')
                            % Re-extract from next point
                            q_R = q(end-2-3*shiftCounter:end-3*shiftCounter);
                            rho_R = q_R(1);
                            en_R = q_R(3);
                            p_R = schm.p(q_R);
                            T_R = p_R / rho_R / Q_;
                            c_R = sqrt(physConst.gamma * Q_ * T_R);
                            M_R = max([0, u_R / c_R]);
                            
                            shiftCounter = shiftCounter + 1;
                            if shiftCounter > 3
                                warning('Shifted a lot to find positive rho, M')
                            end
                        end
                    end
                end

                % Initialize local flags for sonic/subsonic this timestep
                isSonicFlags = [false, false]; % Internal and port
                
                % Capture flow state in PDE domain
                flowState = schm.flowStateR(q);
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    if ~REVERT_MODEL || t < physConst.AirgunCutoffTime
                        if abs(u_R) > 1e-2
                            warning(['Inflow @ t = ' num2str(t) ...
                                '; u|x=0 = ' num2str(u_R)]);
                        end
                    end
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
%                     isSonicFlags(1) = true;
                end
                
                % Compute stagnation pressure
                p0_R = p_R + 0.5 * rho_R * u_R^2;
                
                % reversion: don't trust weak boundary enforcement
                if REVERT_MODEL
                    if M_R >= 1
                        isSonicFlags(1) = true;
                    end
                end

                if REVERT_MODEL
                    if t <= physConst.AirgunCutoffTime
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
                    % Plug flow:
                    % Split velocity at constant pressure to the plug
                    if BYPASS_PLUG
                        uUpstream = u_R;
                    else
                        uUpstream = u_R - shuttle(2);
                    end
%                     
                    
                    
                    %% Case split
                    
                    % Case 0: the port is closed.
                    % Conditions:
                    %   - the port area is zero.
                    % Effects:
                    %   - wall boundary condition applied at the fc exit
                    
                    % Case 1: subsonic everywhere
                    % Conditions:
                    %   - the internal pressure is insufficient to expand
                    %     to the characteristic outflow pressure
                    % Should be an edge case indicating insufficient
                    % internal pressure.
                    % Effects:
                    %   - bubble pressure determines flow state upstream
                    
                    % Case 2: sonic at port
                    % Conditions:
                    %   - Previous case conditions are not met
                    % Should be the most commonly encountered case.
                    
                    % Case 3: sonic at the firing chamber exit
                    % Conditions:
                    %   - the resulting flow to the port must be possible
                    %     (i.e., must be converging in absence of friction)
                    % Effects:
                    %   - the port flow is determined
                    
                    % [Obsolete] Case 4: firing chamber shock formation
                    % Conditions:
                    %   - the flow was internally sonic, but the resulting
                    %     steady flow to the port was not possible
                    % e.g. if the flow is choked inside and the port is
                    % closed too rapidly
                    % Effects:
                    %   - the port becomes sonic
                    %   - shock formation must come through Euler domain BC
                    % Essentially the same as case 4.
                    
                    % Compute critical back pressure at characteristic
                    % bubble velocity

                    % Case 0: port is closed
                    if APortExposed == 0 %|| M_R <= 0
                        caseNum = 0; % Port is closed
                        vel_a = 0;
                        massFlowPort = 0;
                        % Smooth continuation
                        TPort = T_R / temperatureMachFunction(gamma_, 1);
                        pPort = p_R / pressureMachFunction(gamma_, 1);
                        rhoPort = pPort / Q_ / TPort;
%                     elseif p0_R * pressureMachFunction(gamma_, 1) < ...
%                     pBubble % Original issue
                    elseif p0_R < pBubble
                        caseNum = 1; % Subsonic
                        if DEBUG
                            assert(uUpstream >= 0);
                        end          
                        % Mass flow rate is computed using pressure as
                        % influenced by downstream
%                         [pTarget, TPort, rhoPort, cPort, ...
%                             massFlowPort] ...
%                         [pTarget, rhoPort, TPort] ...
%                             = resolveSubsonicPressure(obj, ...
%                                APortExposed, M_R, T_R, pBubble);
                           
                      [MTarget, rhoPort, TPort] = resolveSubsonicVelocity(obj, ...
                               APortExposed, p_R, T_R, pBubble);

                        massFlowPort = rho_R * u_R * ...
                            physConst.crossSectionalArea;
                           
                        % Enforce pressure continuity through port
                        pPort = pBubble;
%                         u_RTarget = vel_a;
                        
                    elseif APortExposed > physConst.crossSectionalArea
                        % Resolve chamber-choked flow
                        caseNum = 3;
                        [velocityPort, pPort, TPort, ...
                            rhoPort, cPort, massFlowPort] ...
                            = resolveSonicChamber(obj, APortExposed, uUpstream/c_R, p_R, T_R);
                        isSonicFlags(1) = true;
                    else
                        % Port-choked flow
                        caseNum = 2;
                        if ~isempty(p_RTarget)
                            p_R_Used = p_RTarget;
                        else
                            p_R_Used = p_R;
                        end
                        [velUp, pPort, TPort, rhoPort, cPort, massFlowPort]...
                        = resolveSonicPort(obj, APortExposed, uUpstream/c_R, p_R_Used, T_R);
                        % Add constraint velocity of plug flow
                        if BYPASS_PLUG
                            vel_a = velUp;
                        else
                            vel_a = velUp + shuttle(2);
                        end
                        u_RTarget = vel_a;
                        isSonicFlags(2) = true;
                    end 
                    
                    % Compute throat pressure
                    if ~isempty(p_RTarget)
                        pThroat = p_RTarget;
                    else
                        pThroat = p_R;
                    end
                    
                    pCrit = 1; % Unused
                    
                    if ~REVERT_MODEL
                        assert(T_R > 0 && ...
                               p_R > 0  && ...
                               c_R > 0);
                    end
                end

                %% Console output
                % Just building a string
                str = sprintf("TIME %.3e", t);
                if isSonicFlags(1)
                    str = str + " | INT    SONIC";
                else
                    str = str + " | INT SUBSONIC";
                end
                if isSonicFlags(2)
                    str = str + " | PORT    SONIC";
                else
                    str = str + " | PORT SUBSONIC";
                end
                % Data
                str = str + " | SHUTPOS: " + sprintf("%.4e", shuttle(1));
                str = str + " | SHUTVEL: " + sprintf("%.4e", shuttle(2));
%                 str = str + " | UPORT*: " + sprintf("%.4e", uPortGuess);
%                 str = str + ...
%                     " | MdotPORT: " + sprintf("%.4e", massFlowPort);
%                 str = str + ...
%                     " | APORT/A*: " + sprintf("%.4e", APortOnASonic);
%                 fprintf(str + "\n");

                %% State regularization
                if M_R < 0 && ~REVERT_MODEL
                    if abs(M_R) > 1e-6
                        warning('Negative mach number');
                    end
                    M_R = 0;
                    u_R = 0;
                end

                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution
                % Early data: should be initial ~80g accel
                if ~REVERT_MODEL
                    volPlug = physConst.plugVolumeFn(shuttle(1));
                    assert(volPlug > 0);
%                     assert(plug(1) > 0);
%                     assert(plug(2) > 0);
                    rhoPlug = plug(1) /  volPlug;
                    TPlug = plug(2) / plug(1) / physConst.c_v;
                    pPlug = rhoPlug * Q_ * TPlug;                           % TODO: output and complete
                    % HACK: override pPlug with p_R
                    pPlug = p_R;
                    % Send boundary pressure to shuttle assembly
                    [dShuttle, pShutRear, pShutFront, pMid, shuttleMonitor] = ...
                        shuttleEvolve(...
                        shuttle, p_R, ...
                        physConst, obj.opChamber);
                else
                    dShuttle = 0*shuttle;
                end
                
                % Check for overshoot at weakly enforced sonic BC
%                 assert(~REVERT_MODEL || u_R < 500)

                %% State monitoring
                if ~REVERT_MODEL
                    assert(all(isreal([M_R, u_R, c_R, ...
                        massFlowPort])));
%                     assert(all([M_R, u_R, c_R, ...
%                         massFlowPort] >= 0));
                    
%                     if APortExposed > 0
%                         uPort = massFlowPort/(rhoPort * APortExposed);
%                     else
%                         uPort = 0;
%                     end
%                     M_port = uPort / sqrt(gamma_ * Q_ * TPort);
%                     monitor(11,1) = p0_R * pressureMachFunction(gamma_, 1) / pBubble;
%                     monitor(12,1) = pCrit;
%                     monitor(14,1) = dShuttle(2) * physConst.shuttleAssemblyMass * shuttle(2);
            

                    if caseNum == 0
                        exportedpTarget = NaN;
                        exporteduTarget = 0;
                    elseif caseNum == 1
                        exportedpTarget = NaN;
                        exporteduTarget = MTarget * c_R;
%                         exportedpTarget = pTarget;
%                         exporteduTarget = NaN;
                    elseif caseNum == 2
                        exportedpTarget = NaN;
                        exporteduTarget = u_RTarget;
                    elseif caseNum == 3
                        exportedpTarget = NaN;
                        exporteduTarget = NaN;
                    else
                        error('Unknown case.');
                    end

                    monitor = struct('t', t, ...
                        'M_R', M_R, ...
                        'M_port', NaN, ... % deprecated
                        'caseNum', caseNum, ...
                        'u_R', u_R, ...
                        'uPort', NaN, ... % deprecated
                        'shuttlePos', shuttle(1), ...
                        'shuttleVel', shuttle(2), ...
                        'shuttleOpRearMass', shuttle(3), ...
                        'shuttleOpRearEnergy', shuttle(4), ...
                        'shuttleOpFrontMass', shuttle(5), ...
                        'shuttleOpFrontEnergy', shuttle(6), ...
                        'pRatio', pThroat/pCrit, ...
                        'ARatio', APortExposed / ...
                                  physConst.crossSectionalArea, ...
                        'p_R', p_R, ...
                        'pShutRear', pShutRear, ...
                        'pShutFront', pShutFront, ...
                        'pMid', pMid, ...
                        'bubbleR', bubble(1), ...
                        'bubbleRDot', bubble(2), ...
                        'bubbleMass', bubble(3), ...
                        'bubbleEnergy', bubble(4), ...
                        'exportedpTarget', exportedpTarget, ...
                        'exporteduTarget', exporteduTarget ...
                        );
                else
                    % Do nothing
                end

                %% Airgun PDE evolution
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                if REVERT_MODEL
                    if velocityPort == 0 && t > 0.0
                        dq = dq + closure_r_closed(q);
%                         assert(t >= 0.010)
                    elseif ~isSonicFlags(1) % Subsonic (in chamber)
                        dq = dq + closure_r_out_sub(q, pBubble);
                    end
                else
                    if caseNum == 0
                        % Closed end when A_port = 0
                        dq = dq + closure_r_closed(q);
                    elseif caseNum == 1
                        % Pressure boundary condition up
%                         dq = dq + closure_r_out_sub(q, pTarget);
                        disp("Subsonic u condition:" + MTarget*c_R);
                        dq = dq + closure_r_out_sub_vel(q, MTarget*c_R);
                    elseif ~isSonicFlags(1) % Subsonic in chamber
                        if vel_a == 0
                            dq = dq + closure_r_closed(q);
                        elseif isSonicFlags(2) % Sonic at port: momentum BC
                            if schm.flowStateR(q) ~= scheme.Euler1d.SUPERSONIC_OUTFLOW
                                dq = dq + closure_r_out_sub_vel(q, vel_a); % TODO: CHECK IF IT'S REAL
                            end
%                             disp(vel_a)
                        else % Subsonic at port: pressure matching
%                             warning('Unchecked subsonic port BC.')
                            dq = dq + closure_r_out_sub_vel(q, vel_a);
%                             dq = dq + closure_r_out_sub(q, pUpstream);
                        end
                    else
                        % Check the sonic claim
%                         assert(schm.flowStateR(q) == scheme.Euler1d.SUPERSONIC_OUTFLOW)
                        %                         warning('Sonic-in-chamber case not treated.')
                        warning('Discr -> RHS logged a sonic case. No BC applied.')
                    end
                end
                
                if ~all(isreal(dq))
                    warning('Complex dq, taking real part')
                    dq = real(dq);
                end
%                 assert(all(isreal(dq)))
                
                %% Bubble evolution
                % Compute port velocity
                if APortExposed > 0 && rhoPort > 0
                    velocityPort = massFlowPort / rhoPort / APortExposed;
                else
                    velocityPort = 0;
                end
                % Compute port specific stagnation energy
                if REVERT_MODEL
                    ePort = en_R;
                else
                    % Compute total energy per volume
                    ePort = rhoPort * physConst.c_v * TPort + ...
                        0.5 * rhoPort * velocityPort^2;
                end
                % Compute bubble differential
                dBubble = bubbleRHS(bubble, ...
                    rhoPort, ...
                    velocityPort, ...
                    ePort, ...
                    pPort, ...
                    APortExposed, ...physConst.A, ...APortExposed, ...
                    physConst);
                % Freeze bubble when port is closed
                if ~REVERT_MODEL
                    if APortExposed == 0 || abs(M_R) < 1e-6
                        if obj.bubbleFrozen
                            dBubble = zeros(size(dBubble));
                        end
                    else
                        obj.bubbleFrozen = false;
                    end
                    
                    assert(dBubble(3) >= 0)
%                     assert(bubble(3) >= 0)
                end
                
                %% Plug evolution
                if REVERT_MODEL
                    dPlug = [0;0];
                else
                    dPlug = [nan; nan];
                    dPlug(1) = rho_R * u_R * physConst.crossSectionalArea;
                    dPlug(2) = dPlug(1) * physConst.c_v * T_R + ...
                        0.5 * dPlug(1) * u_R.^2;
                end
            end
            obj.RHS = @RHS;
        end
        
        
        
        
        % Subsonic port, subsonic firing chamber
        % Bubble parameters matter; information can propagate upstream.
        % Note that at the bubble interface (assumed spherical) the
        % velocity is equal to the bubble expansion velocity \dot{R}.
        % 
        % We assume the expansion is isentropic, and furthermore that
        % shocks are absent (although the subsonic resolution may be
        % invoked right after the float becomes unchoked, and in reality
        % there could still be reflecting shocks present).
        %
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_R            Upstream mach number
        %   p_R            Upstream pressure
        %   T_R            Upstream temperature
        %   ADownstream
        %   MDownstream
        %   pDownstream
        %   TDownstream
        % Output:
        %   vel_a          Upstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        %   pUpstream      Upstream pressure
        % Pre:
        %   APortExposed positive or negative (positive part taken).
        %   M_R > 0
        function [pUpstream, rhoPort, TTotal]..., TPort, rhoPort, cPort, massFlowPort] ...
                 = resolveSubsonicPressure(obj, APortExposed, ...
                MUpstream, TUpstream, pDownstream)
            %% Pre
%             assert(MUpstream >= 0);
            % Workaround: capping the mach number due to overloss to the
            % plug flow
            MUpstream = max([0, MUpstream]);

            %% Take positive part of the exposed area
            APortExposed = max([0, APortExposed]);

            %% Compute upstream pressure from downstream information [V1]
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            
            % Compute the stagnation pressure from output
            pTotal = pDownstream;
            % Compute upstream pressure
            pUpstream = pTotal * pressureMachFunction(gamma_, MUpstream);
            % Compute stagnation energy from upstream
            TTotal = TUpstream / temperatureMachFunction(gamma_, MUpstream);             % Q: Which one to use?
            
            rhoPort = pTotal / (Q_ * TTotal);
            
%             % Compute sonic area from upstream
%             ASonic = obj.physConst.crossSectionalArea / ...
%                 areaMachFunction(gamma_, MUpstream);
%             % Use relative area for M Port
%             MPort = obj.machAreaFunction(APortExposed/ASonic);
%             % Compute the rest of the state at port
%             TPort = TTotal * temperatureMachFunction(gamma_, MPort);
%             cPort = sqrt(gamma_ * Q_ * TPort);
%             pPort = pTotal * pressureMachFunction(gamma_, MUpstream);
%             rhoPort = pPort / (Q_ * TPort);
%             massFlowPort = MPort * cPort * rhoPort * APortExposed;
%             
%             % Compute chamber exit velocity (just unwinding the Mach #)
%             vel_a = MUpstream * sqrt(gamma_ * Q_ * TUpstream);
        end
        
        % !--[Work in progress]--!
        function [MUpstream, rhoPort, TTotal]..., TPort, rhoPort, cPort, massFlowPort] ...
                 = resolveSubsonicVelocity(obj, APortExposed, ...
                pUpstream, TUpstream, pDownstream)
            %% Take positive part of the exposed area
            APortExposed = max([0, APortExposed]);

            %% Compute upstream pressure from downstream information [V1]
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            
            % Compute the stagnation pressure from output
            pTotal = pDownstream;
            % Compute mach number
            MUpstream = machPressureFunction(gamma_, ...
                pUpstream/pTotal);
            
            % Compute stagnation energy from upstream
            TTotal = TUpstream / temperatureMachFunction(gamma_, MUpstream);             % Q: Which one to use?
            
            rhoPort = pTotal / (Q_ * TTotal);
            
%             % Compute sonic area from upstream
%             ASonic = obj.physConst.crossSectionalArea / ...
%                 areaMachFunction(gamma_, MUpstream);
%             % Use relative area for M Port
%             MPort = obj.machAreaFunction(APortExposed/ASonic);
%             % Compute the rest of the state at port
%             TPort = TTotal * temperatureMachFunction(gamma_, MPort);
%             cPort = sqrt(gamma_ * Q_ * TPort);
%             pPort = pTotal * pressureMachFunction(gamma_, MUpstream);
%             rhoPort = pPort / (Q_ * TPort);
%             massFlowPort = MPort * cPort * rhoPort * APortExposed;
%             
%             % Compute chamber exit velocity (just unwinding the Mach #)
%             vel_a = MUpstream * sqrt(gamma_ * Q_ * TUpstream);
        end
        
        
        
        
        % Subsonic port, subsonic firing chamber
        % Bubble parameters matter; information can propagate upstream.
        % Note that at the bubble interface (assumed spherical) the
        % velocity is equal to the bubble expansion velocity \dot{R}.
        % 
        % We assume the expansion is isentropic, and furthermore that
        % shocks are absent (although the subsonic resolution may be
        % invoked right after the float becomes unchoked, and in reality
        % there could still be reflecting shocks present).
        %
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_R            Upstream mach number
        %   p_R            Upstream pressure
        %   T_R            Upstream temperature
        %   ADownstream
        %   MDownstream
        %   pDownstream
        %   TDownstream
        % Output:
        %   vel_a          Upstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        %   pUpstream      Upstream pressure
        % Pre:
        %   APortExposed positive or negative (positive part taken).
        %   M_R > 0
        function [TDownstream, massFlowRate, uUpstream, MUpstream, ...
                  rhoDownstream] ...
                 = resolveSubsonicByArea(obj, AUpstream, p0Upstream, ...
                TUpstream, ADownstream, pDownstream, machAreaFn)
            if ADownstream == 0
                massFlowRate = 0;
                TDownstream = TUpstream;
                uUpstream = 0;
                MUpstream = 0;
                rhoDownstream = pDownstream / ...
                    (Q_ * TDownstream);
                return
            end
            
            % Alias
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            
            % HACK: subpressure
            if pDownstream/p0Upstream >= 1
                error('Downstream pressure is greater thna upstream pressure.')
                p0Upstream = pDownstream;
            end
            
            % Compute downstream M from pressure, p0
            MDownstream = machPressureFunction(gamma_, ...
                pDownstream/p0Upstream);
            ASonic = ADownstream / areaMachFunction(gamma_, MDownstream);
            
            % HACK: subpressure
            if AUpstream/ASonic < 1
                ASonic = AUpstream;
            end
            
            % Use area upstream to determine upstream M
            MUpstream = machAreaFn(AUpstream / ASonic);
            cUpstream = sqrt(gamma_ * Q_ * TUpstream);
            uUpstream = MUpstream * cUpstream;

            % Compute knowing upstream T
            TTotal = TUpstream / ...
                temperatureMachFunction(gamma_, MUpstream);
            
            % Downstream computations
            TDownstream = TTotal * ...
                temperatureMachFunction(gamma_, MDownstream);            
            rhoDownstream = pDownstream / Q_ / TDownstream;
            cDownstream = sqrt(gamma_ * Q_ * TDownstream);
            massFlowRate = rhoDownstream * MDownstream * cDownstream * ...
                ADownstream;
        end
        
        
        
        
        
        
        
        % Sonic port, subsonic firing chamber
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_R            Upstream mach number
        %   p_R            Upstream pressure
        %   T_R            Upstream temperature
        % Output:
        %   vel_a          Upstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        % Pre:
        %   APortExposed positive or negative (positive part taken).
        function [vel_a, pPort, TPort, rhoPort, cPort, massFlowPort]...
                = resolveSonicPort(obj, APortExposed, M_R, p_R, T_R)
            %% Compute velocity boundary condition
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            % Set sonic area to be the positive part of the
            % current port area
            A_sonic = max([0,APortExposed]);            
            
            % Numerical tolerancing
            if obj.physConst.crossSectionalArea / A_sonic > 1 ...
                && obj.physConst.crossSectionalArea / A_sonic < 1.2
                A_sonic = obj.physConst.crossSectionalArea;
            end
            
            % Compute chamber outlet velocity as dependent
            if APortExposed <= 0
                vel_a = 0;
            else
                % Compute upstream mach number (set M_a as
                % boundary condition)
                areaRatio = obj.physConst.crossSectionalArea / A_sonic;
                % Fix the area ratio
                areaRatio = max(areaRatio, 1);
                M_a = obj.machAreaFunction(areaRatio);
                vel_a = M_a * sqrt(gamma_ * Q_ * T_R);
            end            
            %% Compute choked properties at port
            pPort = p_R * pressureMachFunction(gamma_, 1) ...
                / pressureMachFunction(gamma_, M_R);
            TPort = T_R * temperatureMachFunction(gamma_, 1)...
                / temperatureMachFunction(gamma_, M_R);
            rhoPort = pPort / Q_ / TPort;
            cPort = sqrt(gamma_ * Q_ * ...
                TPort);
            massFlowPort = rhoPort * cPort * APortExposed;
            
            %% Sanity check
            assert(TPort > 0 && ...
                pPort > 0 && ...
                rhoPort > 0 && ...
                cPort > 0 && ...
                isreal(cPort))
        end
        
        
        % Subsonic port, sonic firing chamber (subsonic inflow)
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_L            Upstream mach number
        %   p_L            Upstream pressure
        %   T_L            Upstream temperature
        % Output:
        %   velPort        Downstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        % Pre:
        %   APortExposed >= 0
        %   M_L <= 1
        %   p_L > 0
        %   T_L > 0
        function [velPort, pPort, TPort, rhoPort, cPort, massFlowPort] ...
                = resolveSonicChamber(obj, APortExposed, M_L, p_L, T_L)
            %% Precheck
            if M_L > 1
%                 warning('Shuttle velocity acting up!')
                M_L = min(1, M_L);
            end
            assert(APortExposed >= 0 && ...
                   M_L <= 1 && ...
                   p_L > 0 && ...
                   T_L > 0);
            %% Define aliases and constants
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            %% Compute port conditions using forward relations
            ASonic = obj.physConst.crossSectionalArea / ...
                areaMachFunction(gamma_, M_L);
            portAreaRatio = APortExposed/ASonic;
            % Check for shock formation due to limited port size
            if portAreaRatio < 1
                portAreaRatio = 1;
%                 warning('Not really shock formation. Redirect.');
            end
            MPort = obj.machAreaFunction(portAreaRatio);
            pPort = p_L * pressureMachFunction(gamma_, MPort) / ...
                pressureMachFunction(gamma_, M_L);
            TPort = T_L * temperatureMachFunction(gamma_, MPort) / ...
                temperatureMachFunction(gamma_, M_L);
            rhoPort = pPort / Q_ / TPort;
            
            cPort = sqrt(gamma_ * Q_ * TPort);
            velPort = cPort * MPort;            
            massFlowPort = rhoPort * velPort * APortExposed;
        end
    end    
end

