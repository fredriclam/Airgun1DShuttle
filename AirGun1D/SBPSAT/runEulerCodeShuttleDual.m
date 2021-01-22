% Runs coupled Euler codes simultaenously for both the shuttle-coupled
% model and the instant-open, shuttle-free model.
%
% Input:
%   nx -- Number of grid points to use per meter
%   tspan -- Simulation window [s]
%   paramAirgun -- Struct of airgun parameters
%   runShuttleFreeFlag -- (default=false) Whether the shuttle-free model
%                         runs
%   metadata -- (optional) metadata struct to append to
%
% Output:
%   solution -- Complete solution struct for coupled model
%   metadata -- Updated metadata struct
%   solShuttleFree -- Complete solution struct for shuttle-free model.
%                     If runShuttleFreeFlag==false, returns []

function [solution, metadata, solShuttleFree] = ...
    runEulerCodeShuttleDual(nx,tspan,paramAirgun,...
                            runShuttleFreeFlag,metadata)
    %% Defaults for this function
    % Plot every [s]
    PLOT_INTERVAL = 0*1e-4;
    REL_TOL = 1e-4;
                        
    %% User checking and unpacking input parameters
    if nargin < 3 || nargin > 5
        error("Unknown number of parameters provided.");
    elseif nargin == 4
        runShuttleFreeFlag = false;
    end
    
    if ~runShuttleFreeFlag
        % Set empty values for shuttlefree solution
        solShuttleFree = [];
    end
    
    try
        airgunPressure = paramAirgun.airgunPressure;
        airgunLength = paramAirgun.airgunLength;
        airgunCrossSecArea = paramAirgun.airgunCrossSecArea;
        airgunPortArea = paramAirgun.airgunPortArea;
        airgunDepth = paramAirgun.airgunDepth;
        airgunFiringChamberProfile = paramAirgun.airgunFiringChamberProfile;
        airgunOperatingChamberProfile = paramAirgun.airgunOperatingChamberProfile;
        bubbleInitialVolume = paramAirgun.bubbleInitialVolume;
        shuttleBdryPenaltyStrength = paramAirgun.shuttleBdryPenaltyStrength;
    catch
        error("Error unpacking airgun parameters. " ...
              + "Some values may not be specified.");
    end
    
    wallClockStart = tic();

    %% Initialize airgun models (reverted and new model)
    orderSBP = 3; 
    discretization = ...
        DiscrAirgunShuttleMulti(...
            nx, ...
            orderSBP, ...
            airgunPressure, ...
            airgunLength, ...
            airgunCrossSecArea, ...
            airgunDepth, ...
            airgunPortArea, ...
            false, ...
            airgunFiringChamberProfile, ...
            airgunOperatingChamberProfile, ...
            bubbleInitialVolume, ...
            shuttleBdryPenaltyStrength);
    q0 = discretization.q0;
    bubble0 = discretization.bubble0;
    shuttle0 = discretization.shuttle0;
    RHS = discretization.RHS;
    
    % Housekeeping
    metadata.discretization = discretization;
    metadata.q0 = q0;
    metadata.bubble0 = bubble0;
    metadata.shuttle0 = shuttle0;
    
    % Plot counter
    lastPlot = 0;
    
    function dy = odefun(t,y)
        %% Unpack state
        q = y(1:length(q0));
        bubble = y(length(q0)+1:...
                   length(q0)+length(bubble0));
        shuttle = y(length(q0)+length(bubble0)+1:...
                    length(q0)+length(bubble0)+length(shuttle0));
        %% Compute RHS
        [dq, dBubble, dShuttle] = ...
            RHS(q,t,bubble,shuttle);
        %% Repack state
        dy = [dq; dBubble; dShuttle];

        %% Plot every PLOT_INTERVAL
        if t - lastPlot >= PLOT_INTERVAL
            figure(99); clf;
            %% Velocity plot
            subplot(3,1,1);
            u = q(2:3:end) ./ q(1:3:end);
            plot(discretization.schm.x(2:3:end), u);
            ylabel('u [m/s]')
            title(num2str(t))
            hold on
            fs = discretization.fullState(q, t, bubble, shuttle, false);
            if strcmpi(fs.portStates.caseKey, 'portChoked')
                plot(discretization.schm.x(2:3:end), ...
                    fs.portStates.velocityPort*...
                      ones(size(discretization.schm.x(2:3:end))), ...
                    '--r');
            elseif strcmpi(fs.portStates.caseKey, 'chamberChokedForced')
                plot(discretization.schm.x(2:3:end), ...
                    fs.portStates.velocityPort*...
                      ones(size(discretization.schm.x(2:3:end))), ...
                    '--b');
            end
            hold off
            
            %% Temperature plot
            subplot(3,1,2);
            T = (q(3:3:end) - 0.5 * q(2:3:end).^2 ./ q(1:3:end)) ./ ...
                q(1:3:end) ./discretization.physConst.c_v;
            plot(discretization.schm.x(2:3:end), T);
            xlabel('$t$ [s]', 'Interpreter', 'latex')
            ylabel('$T$ [K]', 'Interpreter', 'latex')
            
            %% Pressure plot
            subplot(3,1,3);
            p = (discretization.physConst.gamma - 1) * (q(3:3:end) - 0.5 * q(2:3:end).^2 ./ q(1:3:end));
            plot(discretization.schm.x(2:3:end), p);
            xlabel('$t$ [s]', 'Interpreter', 'latex')
            ylabel('$p$ [Pa]', 'Interpreter', 'latex')
            
            drawnow;
            lastPlot = t;
        end
    end

    y0 = [q0; bubble0; shuttle0];

    % On setup completion
    metadata.wallClockSetupElapsed = toc(wallClockStart);    
    disp('Setup complete. Starting ODE solver.')
    
    %% Run ODE solver
    % Set ODE solver options
    options = odeset('RelTol',REL_TOL);
    sol_ode = ode45(@odefun, tspan, y0,options);

    %% Split quantities from state space vector
    ind1 = 1;
    ind2 = length(q0);
    q = sol_ode.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble = sol_ode.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle = sol_ode.y(ind1:ind2, :);
    
    solution = struct( ...
        'soln', sol_ode, ...
        'q', q, ...
        'bubble', bubble, ...
        'shuttle', shuttle ...
    );
    
    %% Old model [WIP]
    if runShuttleFreeFlag
        dRevert = DiscrAirgunShuttleMulti(nx,orderSBP,airgunPressure,airgunLength, ...
            airgunPortArea,airgunDepth,0,true);
        
        % Grab initial states for each subsystem (same for both systems)
        q0 = dRevert.q0;
        bubble0 = dRevert.bubble0;
        RHS = dRevert.RHS;
        % Grab ODE's RHS objects (Shuttle version has new boundary condition
        % defined for the Euler domain)
        
        if true % in function RHS
            % Extract concatenated data
            ind1 = 1;
            ind2 = ind1+length(q0)-1;
            q2 = y(ind1:ind2);

            ind1 = ind2+1;
            ind2 = ind1+length(bubble0)-1;
            bubble2 = y(ind1:ind2);

            ind1 = ind2+1;
            ind2 = ind1+length(shuttle0)-1;
            shuttle2 = y(ind1:ind2);

            ind1 = ind2+1;
            ind2 = ind1+length(plug0)-1;
            plug2 = y(ind1:ind2);

            % Update new model
            [dq, dBubble, dShuttle, dPlug, p_RTarget, u_RTarget, monitor, shuttleMonitor] = ...
                RHS(q2,t,bubble2,shuttle2,plug2,p_RTarget,u_RTarget);
            dyShuttle = [dq; dBubble; dShuttle; dPlug];
        end
    end
    
    % Use same initial bubble size for both
    y0 = [q0; bubble0; shuttle0];

    %% Run ODE solver
    % Set ODE solver options
    options = odeset('RelTol',1e-4);
    sol_ode = ode45(@odefun, tspan, y0,options);
    
    ind1 = 1;
    ind2 = ind1+length(q0)-1;
    q2 = sol_ode.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble2 = sol_ode.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle2 = sol_ode.y(ind1:ind2, :);
    
    % End wall clock timer
    
    metadata.wallClockTotalElapsed = toc(wallClockStart);
end

