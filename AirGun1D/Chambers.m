%% Chambers object
% Contains data for the middle and operating chamber. 
%% TODO: visualize operating chamber geometry

classdef Chambers
    properties
        totalVolume
        piston_diam
        shuttle_area_right_front
        shuttle_area_right_rear
        total_travel_length
    end
    
    methods
        function obj = Chambers()
            %% Set data for operating chamber
            % Pistom diameter [m]
            obj.piston_diam = 11.1 * 0.0254;
            % Shuttle area in the operating chamber, front [m^2]
            obj.shuttle_area_right_front = pi/4 * ( ...
                obj.piston_diam^2);
            % Shuttle area in the operating chamber, rear [m^2]
            %   Wetted area is smaller due to shaft diameter
            obj.shuttle_area_right_rear = pi/4 * ( ...
                obj.piston_diam^2 - (2.1 * 0.0254)^2); %
            % Total travel length of the shuttle [m]
            obj.total_travel_length = 3.009 * 0.0254;
            
            % Compute constrained geometric parameters
            obj.totalVolume = obj.rearVolume(obj.total_travel_length);
        end
        
        % Visualize the chambers (without the piston position)
        function obj = visualizeStatic(obj)
            % Visualize operating chamber
            xControlPoints = linspace(0, obj.total_travel_length, 100);
            yControlPoints = obj.piston_diam/2 + ...
                             arrayfun(@obj.gapProfile, xControlPoints);
            
            loopPoints = @(z) [z, fliplr(z)];
            mirrorXAxis = @(z) [-z, fliplr(z)];
            
            xMirrored = loopPoints(xControlPoints);
            yMirrored = mirrorXAxis(yControlPoints);
            
            
            fill(xMirrored, yMirrored,  [200 200 200]/255)
            axis equal
        end
        
        % Visualize the chambers with piston position
        function obj = visualize(obj, xi)
            error('WIP');
            obj.visualizeStatic();
        end
         
        % Return the volume of gas in front of the piston
        function V = frontVolume(obj, x)
            V = obj.totalVolume - obj.rearVolume(x);
        end
        
        % Return distance between the piston (fixed diameter)
        % and the chamber wall as a function of shuttle displacement [m].
        function gapLength = gapProfile(obj, xi)
            x = xi/obj.total_travel_length;
            if x == 0
                y = 0;
            elseif x < 24/58
                y = 5/58 * x/(24/58);
            elseif x < 35/58
                y = 5/58;
            elseif x < 46/58
                y = 5/58 - 5/58*(x-35/58)/(46/58-35/58);
            elseif x <= 1
                y = 0;
            elseif x < 0 || x > 1
                y = 0;
                warning('Queried geometry outside of bounds')
            end
            gapLength = y * obj.total_travel_length;
        end
        
        % Return area of annular gap, between the piston (fixed diameter)
        % and the chamber wall (varies with shuttle displacement [m]).
        function A = gapArea(obj, xi)
            rOuter = obj.piston_diam/2 + obj.gapProfile(xi);
            rInner = obj.piston_diam/2;
            A = pi*(rOuter^2 - rInner^2);
        end
        
        % Return volume of gas behind the shuttle, along the
        % shuttle travel axis, as a function of the shuttle
        % displacement [m]. 
        function V = rearVolume(obj, xi)
            % Returns rear volume as function of xi
            x = xi/obj.total_travel_length;
            % Analytic integral of gap profile as function of x/l
            % Using formula for conical frustum
            if x == 0
                V = 0;
            elseif x < 24/58
                % Conic frustrum
                frustumRadius1 = obj.piston_diam/2;
                frustumRadius2 = obj.piston_diam/2 + obj.gapProfile(xi);
                V = pi/3 * xi * (frustumRadius1^2 + ...
                    frustumRadius1*frustumRadius2 + ...
                    frustumRadius2^2);
            elseif x < 35/58
                % Conic frustrum plus...
                frustumRadius1 = obj.piston_diam/2;
                frustumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustumRadius1^2 + ...
                    frustumRadius1*frustumRadius2 + ...
                    frustumRadius2^2);
                % ...washer
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (x-24/58) * obj.total_travel_length;
            elseif x < 46/58
                % Conic frustrum plus...
                frustumRadius1 = obj.piston_diam/2;
                frustumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustumRadius1^2 + ...
                    frustumRadius1*frustumRadius2 + ...
                    frustumRadius2^2);
                % ...washer plus...
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (35/58-24/58) * obj.total_travel_length;
                % ...conic frustrum
                frustumRadius1 = obj.piston_diam/2 + ...
                    obj.gapProfile(35/58 * obj.total_travel_length);
                frustumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(xi);
                V = V + pi/3 * (x-35/58) * obj.total_travel_length ...
                    * (frustumRadius1^2 + ...
                    frustumRadius1 * frustumRadius2 + ...
                    frustumRadius2^2);
            elseif x <= 1
                % Conic frustrum plus...
                frustumRadius1 = obj.piston_diam/2;
                frustumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustumRadius1^2 + ...
                    frustumRadius1*frustumRadius2 + ...
                    frustumRadius2^2);
                % ...washer plus...
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (35/58-24/58) * obj.total_travel_length;
                % ...conic frustrum plus...
                frustumRadius1 = obj.piston_diam/2 + ...
                    obj.gapProfile(35/58 * obj.total_travel_length);
                frustumRadius2 = obj.piston_diam/2;
                V = V + pi/3 * (46/58-35/58) * obj.total_travel_length ...
                    * (frustumRadius1^2 + ...
                    frustumRadius1 * frustumRadius2 + ...
                    frustumRadius2^2);
                % ...washer
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(58/58*obj.total_travel_length))^2) * ...
                    (x-46/58) * obj.total_travel_length;
                
            elseif x < 0 || x > 1
                V = 0;
                warning('Queried geometry outside of bounds')
            end
        end
        
        % Generate static plot of volume behind the shuttle, along the
        % shuttle travel axis.
        function plotVolumeFn(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.rearVolume(xVec(i));
            end
            plot(xVec, yVec, '-', 'LineWidth', 1.5);
            xlabel 'Axis [m]'
            ylabel 'Volume [m^3]'
        end
        
        % Generate static plot of gap length along the shuttle travel axis.
        % Gap allows the gas to flow between the front and the rear of
        % the piston.
        function plotGapLength(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.gapProfile(xVec(i));
            end
             plot(xVec, yVec, '-', 'LineWidth', 1.5);
            xlabel 'Axis [m]'
            ylabel 'Gap length [m]'
        end
        
        % Generate static plot of gap area along the shuttle travel axis.
        % Gap area allows the gas to flow between the front and the rear of
        % the piston.
        function plotGapArea(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.gapArea(xVec(i));
            end
            plot(xVec, yVec, '-', 'LineWidth', 1.5);
            xlabel 'Axis [m]'
            ylabel 'Gap area [m^2]'
        end
    end
end