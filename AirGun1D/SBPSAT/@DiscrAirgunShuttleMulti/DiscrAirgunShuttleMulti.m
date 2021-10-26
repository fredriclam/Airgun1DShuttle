classdef DiscrAirgunShuttleMulti < DiscrAirgun
    properties
        shuttle0              % Initial shuttle state [pos; vel]
        plug0                 % Initial plug control volume state [mass; en]
        machAreaFunction      % Precomputed M(A/A*) function interpolant
        chambers              % Object for middle and operating chambers
        bubbleFrozen          % Freeze bubble state until first open to prevent negative bubble pressure
        wallClockSetupElapsed % Setup timer
        bubbleModel
    end
        methods
        function obj = DiscrAirgunShuttleMulti(nx,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth, airgunCrossSectionalArea, REVERT_MODEL, ...
                airgunFiringChamberProfile, ...
                airgunOperatingChamberProfile, bubbleInitialVolume, ...
                shuttleBdryPenaltyStrength, midChamberMode, ...
                airgunPortLength, bubbleModel, extraOptions)
            
            %% Perform user checks
            if nargin == 8 && ~REVERT_MODEL
                error('Incorrect # of arguments for shuttle-included model.')
            end
            
            %% Default parameters for this function
            DEBUG = false;
            % Whether to use a fixed area test problem instead
            FIXED_AREA_TEST_PROBLEM = false;
            
            %% Default airgun cutoff
            airgunCutoffTime = 0.100;
            
            %% Initial setup
            % Start Timer
            wallClockSetupStart = tic();
            % Call parent constructor
            obj = obj@DiscrAirgun(nx,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            % Add mid and operating chambers
            obj.chambers = Chambers(midChamberMode);
            
            % Add additional configuration with backward-compatible setting
            [physConst, ~, icAirgun, icBubble] = ...
                configAirgun(...
                    'GeneralAirgun', ...
                    airgunPressure, ...
                    airgunLength, ...
                    airgunPortArea, ...
                    airgunDepth, ...
                    airgunCrossSectionalArea, ...
                    airgunFiringChamberProfile, ...
                    airgunOperatingChamberProfile, ...
                    bubbleInitialVolume, ...
                    shuttleBdryPenaltyStrength, ...
                    airgunPortLength, ...
                    extraOptions);
            obj.physConst = physConst;
            obj.physConst.airgunCutoffTime = airgunCutoffTime;
            obj.bubbleModel = bubbleModel;
            if ~REVERT_MODEL
                % Replace this object's description
                obj.description = ...
                    'Airgun augmented with shuttle + mid and op chambers.';
            else
                % Replace this object's description
                obj.description = ...
                    'Airgun with inactive shuttle + mid and op chambers.';
            end
            
            % Alias commonly used objects
            physConst = obj.physConst;
            schm = obj.schm;
            gamma_ = physConst.gamma;
            
            %% Set initial states
            % Replace state from Discr constructor with config values
            obj.q0(1:3:end) = icAirgun.rho0(obj.schm.u);
            obj.q0(2:3:end) = icAirgun.rv0(obj.schm.u);
            obj.q0(3:3:end) = icAirgun.e0(obj.schm.u);
            
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rhov0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);     % Volumetric stagnation energy e
            T0 = (e0-0.5*rhov0^2/rho0)/physConst.c_v/ rho0; % Temperature
            
            % Define initial shuttle state: position; velocity
            % NOTE: position non-zero or else singular--for CV treatment
            if REVERT_MODEL
                xi_0 = 0;
                gasMassPartitionRatio = 1;
                massOpChamber = 0;
            else
                xi_0 = 0.1e-3;
                gasMassPartitionRatio = ...
                    obj.chambers.rearVolume(xi_0) / ...
                    obj.chambers.totalVolume();
                massOpChamber = rho0*obj.chambers.totalVolume();
            end
            
            % Set shuttle initially at rest
            xi_vel_0 = 0;
            massRear = gasMassPartitionRatio * massOpChamber;
            massFront = (1-gasMassPartitionRatio) * massOpChamber;
            volMidChamber0 = physConst.midChamberLength * ...
                physConst.midChamberArea;
            massMid = (physConst.p_inf) / (physConst.Q * T0) * ...
                volMidChamber0;
            
            % Set initial shuttle state vector
            obj.shuttle0 = [
                xi_0; % [m] -- must give rear chamber some room
                xi_vel_0;
                massRear;
                physConst.c_v*T0*massRear;
                massFront;
                physConst.c_v*T0*massFront;
                massMid;
                physConst.c_v*T0*massMid;
            ];   % [m/s]
        
            %% Bubble
            % Freeze bubble until the port first opens (xi >= xi_t)
            obj.bubbleFrozen = true;
            % Replace bubble using provided parameters
            obj.bubble0 = [
                icBubble.R;
                icBubble.Rdot;
                icBubble.m;
                icBubble.E;
            ];

            %% Create boundary condition operators
            % Wall BC
            closure_l = schm.boundary_condition('l', 'wall');
            % Outflow with pressure (for unchoked-everywhere flow)
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');
            % Outflow with velocity (for choked flow)
            closure_r_out_sub_vel = ...
                schm.boundary_condition('r', 'outflow_vel');
            % Wall on the right
            closure_r_closed = schm.boundary_condition('r', 'wall');
            % Characteristic boundary on the right
            closure_r_char = schm.boundary_condition('r', 'char');
            
            closure_r_out_rho = schm.boundary_condition('r', 'outflow_rho');
            
            %% Precompute
            % Precompute mach area function M(A/A*) and, for subsonic,
            % the mach pressure function M(p/p0)
            if ~REVERT_MODEL
                obj.machAreaFunction = precomputeMachAreaFunction(gamma_);
            end
            
            %% Housekeeping
            obj.wallClockSetupElapsed = toc(wallClockSetupStart);
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle] = ...
                    RHS(q, t, bubble, shuttle)
                % RHS function for ode solver
                
                %% Compute state dependent variables
                agState = obj.fullState(q, t, bubble, shuttle, REVERT_MODEL);

                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution
                % Early data: should be initial ~80g accel
                if ~REVERT_MODEL
                    % TODO: fix interface for shuttle
                    % Send boundary pressure to shuttle assembly
                    [dShuttle, pShutRear, pShutFront, pMid, ~] = ...
                        shuttleEvolve(...
                        ...shuttle, agState.eulerDomainStates.p_R, ...
                        shuttle, agState.portStates.pPort, ...
                        physConst, obj.chambers);
                else
                    dShuttle = 0*shuttle;
                end

                %% Set Airgun boundary value problem evolution
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                
                % Get correct right BC closure
                if REVERT_MODEL
                    if t > airgunCutoffTime
                        dq = dq + closure_r_closed(q);
                    elseif schm.flowStateR(q) == scheme.Euler1d.SUBSONIC_INFLOW
                        dq = 0*dq;
                    elseif schm.flowStateR(q) == scheme.Euler1d.SUPERSONIC_OUTFLOW
                        dq = dq; %#ok<ASGSL>
                    else
                        dq = dq + closure_r_out_sub(q, agState.bubbleStates.p);
                    end
                elseif FIXED_AREA_TEST_PROBLEM
                    % Test case: fixed port area override, such that the
                    % flow should be port-choked only
                    if t < 0
                        dq = dq + closure_r_closed(q);
                    else
                        agState = fullState(...
                            obj, q, t, bubble, [1.1*obj.physConst.portLead; 0], REVERT_MODEL);
                        dq = dq + closure_r_out_sub_vel(q, agState.portStates.velocityPort);
                    end
                else
                    % Working code
                    if strcmpi('portClosed', ...
                            agState.portStates.caseKey)
                        dq = dq + closure_r_closed(q);
                    elseif strcmpi('subsonic', ...
                            agState.portStates.caseKey)
                        dq = dq + closure_r_out_sub(q, ...
                            agState.portStates.pPort);
                    elseif strcmpi('portChoked', ...
                            agState.portStates.caseKey)
                        if agState.portStates.velocityPort <= 0
                            % Safeguard backward flow
                            warning('port choked back flow')
                            dq = dq + closure_r_closed(q);
                        else
                            dq = dq + closure_r_out_sub(q, agState.portStates.pPort);
                        end
                    elseif strcmpi('chamberChokedForced', ...
                        agState.portStates.caseKey)
                        % Makes sure transition from port choked to chamber
                        % choked leaves the chamber choked
                        dq = dq + closure_r_out_sub(q, agState.portStates.pPort);
                    end
                end
                
                %% Bubble evolution
                % Compute bubble differential
                
                if REVERT_MODEL
                    dBubble = bubbleRHS(t, bubble, ...
                        agState.portStates.rhoPort, ...
                        agState.portStates.velocityPort, ...
                        agState.portStates.eTotalPort, ...
                        agState.portStates.pPort, ...
                        ... min(airgunCrossSectionalArea*0.0254^2, ...
                        ...    airgunPortArea*0.0254^2), ...
                        airgunCrossSectionalArea*(0.0254)^2, ...
                        physConst, ...
                        bubbleModel ...
                    );
                else
                    dBubble = bubbleRHS(t, bubble, ...
                        agState.portStates.rhoPort, ...
                        agState.portStates.velocityPort, ...
                        agState.portStates.eTotalPort, ...
                        agState.portStates.pPort, ...
                        airgunCrossSectionalArea*(0.0254)^2, ...
                        physConst, ...
                        bubbleModel ...
                    );
                end
                % Freeze bubble when port is closed
                if ~REVERT_MODEL
                    if agState.portStates.APortExposed == 0 ...
                            || abs(agState.eulerDomainStates.M_R) < 1e-6 % TODO: remove M_R condition
                        if obj.bubbleFrozen
                            dBubble = zeros(size(dBubble));
                        end
                    else
                        obj.bubbleFrozen = false;
                    end
                end
            end
            
            % Export RHS function
            obj.RHS = @RHS;
        end
        
        % Declare state computation function
        agState = fullState(obj, q, t, bubble, shuttle, REVERT_MODEL, INCLUDE_ALL_PRIMITIVES)
        
        % Declare constraint enforcement function
        [qTarget, exitFlag] = enforceScalarConstraint(obj, essentialConstraint, q_R)
        
        function pStar = sonicPressure(obj, pressure, M)
            % Compute sonic pressure at given state
            gamma_ = obj.physConst.gamma;
            pStar = pressure * ...
                (((gamma_ + 1)/2) * ...
                (1 + (gamma-1)/2*M.^2)).^(-gamma/(gamma-1));
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
        
        function [MUpstream, exception] = mapChokedPortAreaRatioToM(...
                obj, APortExposed)
            exception = [];
            
            % Take exposed port area as the sonic area
            A_sonic = APortExposed;
            areaRatio = obj.physConst.crossSectionalArea / A_sonic;
            if areaRatio < 1
                exception = struct(...
                    'message', ...
                    'Critical: area ratio is invalid for this case.', ...
                    'details', ['Area ratio must be >= 1 assuming' ...
                                'subsonic upstream (contraction to' ...
                                'sonic)'], ...
                    'areaRatio', areaRatio ...
                );
                MUpstream = 0;
                % Return soft exception and snap to area ratio = 1
                % This allows temporary violation of the relation
                areaRatio = 1;
                MUpstream = obj.machAreaFunction(areaRatio);
            else
                MUpstream = obj.machAreaFunction(areaRatio);
            end
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
            
%             % Numerical tolerancing
%             if obj.physConst.crossSectionalArea / A_sonic > 1 ...
%                 && obj.physConst.crossSectionalArea / A_sonic < 1.2
%                 A_sonic = obj.physConst.crossSectionalArea;
%             end
            
            % Compute chamber outlet velocity as dependent
            if APortExposed <= 0
                vel_a = 0;
            else
                % Compute upstream mach number (set M_a as
                % boundary condition)
                areaRatio = obj.physConst.crossSectionalArea / A_sonic;
%                 % Fix the area ratio
%                 areaRatio = max(areaRatio, 1);
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

