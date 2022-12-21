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

function [solution, metadata] = ...
    runEulerCodeShuttleDual(nx,tspan,paramAirgun,...
                            coupleToShuttle,metadata)
    %% Defaults for this function
    % Plot every [s]
    PLOT_INTERVAL = 1e-2;
    REL_TOL = 1e-5; % 1e-3 to 1e-5 OK
    
    % Final t to run the (uncoupled) bubble model to
    tBubbleFinal = 10;
                        
    %% User checking and unpacking input parameters
    if nargin < 3 || nargin > 5
        error("Unknown number of parameters provided.");
    elseif nargin == 4
        % Default: use coupled shuttle
        coupleToShuttle = true;
    end
    
    try
        airgunPressure = paramAirgun.airgunPressure;
        airgunLength = paramAirgun.airgunLength;
        airgunCrossSecArea = paramAirgun.airgunCrossSecAreaSqInch;
        airgunPortArea = paramAirgun.airgunPortAreaSqInch;
        airgunDepth = paramAirgun.airgunDepth;
        airgunFiringChamberProfile = paramAirgun.airgunFiringChamberProfile;
        airgunOperatingChamberProfile = paramAirgun.airgunOperatingChamberProfile;
        bubbleInitialVolume = paramAirgun.bubbleInitialVolume;
        shuttleBdryPenaltyStrength = paramAirgun.shuttleBdryPenaltyStrength;
        midChamberMode = paramAirgun.midChamberMode;
        airgunPortLength = paramAirgun.airgunPortLength;
        bubbleModel = paramAirgun.bubbleModel;
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
            airgunPortArea, ...
            airgunDepth, ...
            airgunCrossSecArea, ...
            ~coupleToShuttle, ...
            airgunFiringChamberProfile, ...
            airgunOperatingChamberProfile, ...
            bubbleInitialVolume, ...
            shuttleBdryPenaltyStrength, ...
            midChamberMode, ...
            airgunPortLength, ...
            bubbleModel, ...
            metadata.extraOptions);
    q0 = discretization.q0;
    bubble0 = discretization.bubble0;
%     % Lifting to V + V_body
%     V_body = 0.3^3 * pi/ 4; % Approximate extra volume of body included
%     V_new = 4/3*pi*bubble0(1).^3 + V_body;
%     bubble0(1) = (V_new / (4/3*pi))^(1/3);
    
    % Injection for partitioned-energy bubble model
    if isstr(bubbleModel)
        if strcmpi('partition', bubbleModel)
            bubble0 = [bubble0; 0];
        elseif strcmpi('single-power', bubbleModel)
            bubble0(3:4) = bubble0(3:4) / 0.25;
            bubble0 = [bubble0; 0];
        end
    elseif isstruct(bubbleModel)
        if strcmpi('partition', bubbleModel.type)
            bubble0 = [bubble0; 0];
        elseif strcmpi('single-power', bubbleModel.type)
            bubble0(3:4) = bubble0(3:4) / 0.25;
            bubble0 = [bubble0; 0];
        elseif isfield(bubbleModel, 'waterProperties')
            bubble0 = [bubble0; 0];
        end 
    else
%         bubble0 = [bubble0; 0];
        error('Unknown bubble model')
    end
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
        if coupleToShuttle
            shuttle = y(length(q0)+length(bubble0)+1:...
                        length(q0)+length(bubble0)+length(shuttle0));
        else
            shuttle = zeros(size(length(shuttle0)));
        end
        
        %% Compute RHS
        [dq, dBubble, dShuttle] = ...
            RHS(q,t,bubble,shuttle);
        
        %% Repack state
        if coupleToShuttle
            dy = [(t>=0)*dq; dBubble; (t>=0)*dShuttle];
        else
            dy = [(t>=0)*dq; dBubble];
        end

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

    if coupleToShuttle
        y0 = [q0; bubble0; shuttle0];
    else
        y0 = [q0; bubble0];
    end

    % On setup completion
    metadata.wallClockSetupSeconds = toc(wallClockStart);    
    disp('Setup complete. Starting ODE solver.')
    
    %% Run ODE solver
    % Set ODE solver options
    options = odeset('RelTol',REL_TOL);
    figure(99);
    % Compute max wave speed
    gamma = discretization.physConst.gamma;
    h = discretization.h;
    w_max0 = 2 * max( ...
        sqrt(gamma * (gamma - 1) * ( ...
        q0(3:3:end) ./ q0(1:3:end) ...
        - 0.5*q0(2:3:end).^2 ./ q0(1:3:end).^2)));
    % Split ode solve
    if tspan(1) < 0
        % ODE splitting to prevent stiffness in dy
        tspan1 = [tspan(1), 0];
        tspan2 = [0, tspan(2)];
        sol_ode = ode45(@odefun, tspan1, y0,options);
        options = odeset('RelTol',REL_TOL,'MaxStep',h/w_max0);
        sol_ode = odextend(sol_ode, @odefun, tspan2(2), sol_ode.y(:,end), options);
    else
        options = odeset('RelTol',REL_TOL,'MaxStep',h/w_max0);
        sol_ode = ode45(@odefun, tspan, y0,options);
    end

    fprintf('ODE solve complete for coupled phase ([%.3f, %.3f] s).\n', ...
        tspan(1), tspan(2));
    
    %% Split quantities from state space vector
    ind1 = 1;
    ind2 = length(q0);
    q = sol_ode.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble = sol_ode.y(ind1:ind2, :);

    if coupleToShuttle
        ind1 = ind2+1;
        ind2 = ind1+length(shuttle0)-1;
        shuttle = sol_ode.y(ind1:ind2, :);
    end
    
    clear ind1 ind2;
    
    %% Uncoupled bubble solve
    % Provided the provided tspan allows for solution of the coupled system
    % up to the port closing permanently, the following solve takes into
    % account the bubble with no mass or energy input.
    % Reduced memory and computational cost.
    disp('Starting ODE solve for uncoupled bubble.');
    
    solnBubbleContinuation = ode45(...
        @(t,y) bubbleRHS(t, y, 0, 0, 0, 0, 0, discretization.physConst, ...
        paramAirgun.bubbleModel), ...
        [tspan(2), tBubbleFinal], bubble(:,end),options);
    
    bubbleContinuationTime = [sol_ode.x, ...
        solnBubbleContinuation.x(2:end)];
    bubbleContinuationState = [sol_ode.y(length(q0)+1:length(q0)+length(bubble0), :), ...
        solnBubbleContinuation.y(:,2:end)];
    
    %% Data packing
    solution = struct( ...
        'soln', sol_ode, ...
        'q', q, ...
        'bubble', bubble, ...
        'shuttle', shuttle, ...
        'bubbleContinuationTime', bubbleContinuationTime, ...
        'bubbleContinuationState', bubbleContinuationState, ...
        'solnBubbleContinuation', solnBubbleContinuation ...
    );
    
    % End wall clock timer
    metadata.wallClockTotalSeconds = toc(wallClockStart);
    disp('Job finished.')
end

