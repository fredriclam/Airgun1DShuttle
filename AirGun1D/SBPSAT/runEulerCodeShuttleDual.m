% Runs coupled Euler codes simultaenously for both the shuttle-coupled
% model and the instant-open, shuttle-free model.

function [sol, q1, bubble1, shuttle1, plug1, ...
    q2, bubble2, shuttle2, plug2, monitorStates] = ...
    runEulerCodeShuttleDual(nx,airgunPressure,...
    airgunLength,airgunCrossSecArea,airgunPortArea,airgunDepth, ...
    airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength)

    % Initialize airgun models (reverted and new model)
    orderSBP = 3;
    dRevert = DiscrAirgunShuttleMulti(nx,orderSBP,airgunPressure,airgunLength, ...
        airgunPortArea,airgunDepth,0,true);
    dNew = DiscrAirgunShuttleMulti(nx,orderSBP,airgunPressure,airgunLength, ...
        airgunCrossSecArea,airgunDepth,airgunPortArea,false, ...
        airgunFiringChamberProfile, ...
        airgunOperatingChamberProfile, bubbleInitialVolume, ...
        shuttleBdryPenaltyStrength);
    
    % Grab initial states for each subsystem (same for both systems)
    q0 = dRevert.q0;
    bubble0 = dRevert.bubble0;
    shuttle0 = dNew.shuttle0;
    plug0 = dNew.plug0;
    % Grab ODE's RHS objects (Shuttle version has new boundary condition
    % defined for the Euler domain)
    RHSRevert = dRevert.RHS;
    RHSShuttle = dNew.RHS;
    
    % Plot counter
    lastPlot = 0;
    % States to monitor
    clear monitorStates;% = [];
    % Target pressure and u for shuttle-included model
    p_RTarget = [];
    u_RTarget = [];
    function dy = odefun(t,y)
        % Partition y into each subsystem
        ind1 = 1;
        ind2 = length(q0);
        q1 = y(ind1:ind2);

        ind1 = ind2+1;
        ind2 = ind1+length(bubble0)-1;
        bubble1 = y(ind1:ind2);

        ind1 = ind2+1;
        ind2 = ind1+length(shuttle0)-1;
        shuttle1 = y(ind1:ind2);

        ind1 = ind2+1;
        ind2 = ind1+length(plug0)-1;
        plug1 = y(ind1:ind2);
        
        % Update reverted model
        [dq, dBubble, dShuttle, dPlug, ~] = ...
            RHSRevert(q1,t,bubble1,shuttle1,plug1);
        dyRevert = [dq; dBubble; dShuttle; dPlug];

        % Extract concatenated data
        ind1 = ind2+1;
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
        [dq, dBubble, dShuttle, dPlug, p_RTarget, u_RTarget, monitor] = ...
            RHSShuttle(q2,t,bubble2,shuttle2,plug2,p_RTarget,u_RTarget);
        dyShuttle = [dq; dBubble; dShuttle; dPlug];
        
        % Concatenate dy
        dy = [dyRevert; dyShuttle];
        
        % DEBUG: state monitoring
        if ~exist('monitorStates', 'var')
            monitorStates = struct('length', 1000, ...
                'index', 2);
            monitorStates.data = monitor;
            % Duplicate with filler data
            monitorStates.data(1:monitorStates.length) = monitor;
        else
            if monitorStates.index > monitorStates.length
                monitorStates.data(monitorStates.length+1: ...
                                   monitorStates.length+1000) = monitor;
                monitorStates.length = monitorStates.length + 1000;
            end
            monitorStates.data(monitorStates.index) = monitor;
            monitorStates.index = monitorStates.index + 1;
        end

        %% Piggyback plotting every millisecond
        if t - lastPlot >= 1e-3
            figure(99); clf;
            subplot(3,1,1);
            v1 = q1(2:3:end) ./ q1(1:3:end);
            plot(dRevert.schm.x(2:3:end), v1);
            hold on
            v2 = q2(2:3:end) ./ q2(1:3:end);
            plot(dNew.schm.x(2:3:end), v2);
            
            ylabel('u [m/s]')
            title(num2str(t))
            
            subplot(3,1,2);
            T1 = (q1(3:3:end) - 0.5 * q1(2:3:end).^2 ./ q1(1:3:end)) ./ q1(1:3:end) ./dRevert.physConst.c_v;
            plot(dRevert.schm.x(2:3:end), T1);
            hold on
            T2 = (q2(3:3:end) - 0.5 * q2(2:3:end).^2 ./ q2(1:3:end)) ./ q2(1:3:end) ./dNew.physConst.c_v;
            plot(dNew.schm.x(2:3:end), T2);
            xlabel('$t$ [s]', 'Interpreter', 'latex')
            ylabel('$T$ [K]', 'Interpreter', 'latex')
            legend({'Const $A_{port}$', 'Full model'}, 'Interpreter', 'latex', ...
                    'Location', 'eastoutside')
            
            subplot(3,1,3);
            p1 = (dRevert.physConst.gamma - 1) * (q1(3:3:end) - 0.5 * q1(2:3:end).^2 ./ q1(1:3:end));
            plot(dRevert.schm.x(2:3:end), p1);
            ylabel('p [Pa]')
%             equilibriumPressure = ...
%                 dRevert.physConst.shuttle_area_right / ...
%                 dRevert.physConst.cross_sectional_area * ...
%                 dRevert.physConst.p_R0;
            hold on
%             plot(dRevert.schm.x(2:3:end), equilibriumPressure*ones(size(dRevert.schm.x(2:3:end))));
            
            p2 = (dNew.physConst.gamma - 1) * (q2(3:3:end) - 0.5 * q2(2:3:end).^2 ./ q2(1:3:end));
            plot(dNew.schm.x(2:3:end), p2);
            
            drawnow;
            lastPlot = t;
        end
    end

    % Use same initial bubble size for both
    bubble0revert = bubble0;
    y0 = [q0; bubble0; shuttle0; plug0; q0; bubble0revert; shuttle0; plug0];
    
    % Simulation [tmin, tmax]
    % Used values:
    % 10, 30, [150], [210], 300, 400, 600 ms
    tspan = [0; 0.100]; 

    % Set ODE solver options
%     options = odeset('RelTol',1e-3);%, 'MaxStep', 1e-3);
%     options = odeset('RelTol',1e-2);
    options = odeset('RelTol',1e-4);
    
    %% ODE toolbox solver
    sol = ode45(@odefun, tspan, y0,options);
    %% Forward Euler stepper
%     dt = 0.000002;
%     fe_time = tspan(1):dt:tspan(end);
%     y = y0;
%     yHistory = zeros(length(y0), length(fe_time));
%     yHistory(:,1) = y0;
%     for i = 2:length(fe_time)
%         t = fe_time(i);
%         y = y + dt * odefun(t,y);
%         yHistory(:,i) = y;
%     end
%     sol = struct('y', yHistory);

    %% Split quantities from state space vector
    ind1 = 1;
    ind2 = length(q0);
    q1 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble1 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle1 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(plug0)-1;
    plug1 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(q0)-1;
    q2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle2 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(plug0)-1;
    plug2 = sol.y(ind1:ind2, :);
end
