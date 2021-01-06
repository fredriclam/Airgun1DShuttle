classdef OperatingChamber
    properties
        totalVolume
        piston_diam
        shuttle_area_right_front
        shuttle_area_right_rear
        total_travel_length
    end
    
    methods
        function obj = OperatingChamber()
            obj.piston_diam = 11.1 * 0.0254; % [m]
            obj.shuttle_area_right_front = pi/4 * ( ...
                obj.piston_diam^2); % [m^2] 
            obj.shuttle_area_right_rear = pi/4 * ( ...
                obj.piston_diam^2 - (2.1 * 0.0254)^2); % [m^2]
            obj.total_travel_length = 3.009 *0.0254;
            obj.totalVolume = obj.rearVolume(obj.total_travel_length);
        end
         
        function V = frontVolume(obj, x)
            V = obj.totalVolume - obj.rearVolume(x);
        end
        
        function gapLength = gapProfile(obj, xi)
            % Returns deviation from piston edge as function of xi
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
        
        function A = gapArea(obj, xi)
            % Returns total area of annular gap
            rOuter = obj.piston_diam/2 + obj.gapProfile(xi);
            rInner = obj.piston_diam/2;
            A = pi*(rOuter^2 - rInner^2);
        end
        
        function V = rearVolume(obj, xi)
            % Returns rear volume as function of xi
            x = xi/obj.total_travel_length;
            % Analytic integral of gap profile as function of x/l
            % Using formulla for conical frustum
            if x == 0
                V = 0;
            elseif x < 24/58
                % Conic frustrum
                frustrumRadius1 = obj.piston_diam/2;
                frustrumRadius2 = obj.piston_diam/2 + obj.gapProfile(xi);
                V = pi/3 * xi * (frustrumRadius1^2 + ...
                    frustrumRadius1*frustrumRadius2 + ...
                    frustrumRadius2^2);
            elseif x < 35/58
                % Conic frustrum plus...
                frustrumRadius1 = obj.piston_diam/2;
                frustrumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustrumRadius1^2 + ...
                    frustrumRadius1*frustrumRadius2 + ...
                    frustrumRadius2^2);
                % ...washer
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (x-24/58) * obj.total_travel_length;
            elseif x < 46/58
                % Conic frustrum plus...
                frustrumRadius1 = obj.piston_diam/2;
                frustrumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustrumRadius1^2 + ...
                    frustrumRadius1*frustrumRadius2 + ...
                    frustrumRadius2^2);
                % ...washer plus...
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (35/58-24/58) * obj.total_travel_length;
                % ...conic frustrum
                frustrumRadius1 = obj.piston_diam/2 + ...
                    obj.gapProfile(35/58 * obj.total_travel_length);
                frustrumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(xi);
                V = V + pi/3 * (x-35/58) * obj.total_travel_length ...
                    * (frustrumRadius1^2 + ...
                    frustrumRadius1 * frustrumRadius2 + ...
                    frustrumRadius2^2);
            elseif x <= 1
                % Conic frustrum plus...
                frustrumRadius1 = obj.piston_diam/2;
                frustrumRadius2 = obj.piston_diam/2 + ...
                    obj.gapProfile(24/58 * obj.total_travel_length);
                V = pi/3 * 24/58 * obj.total_travel_length * (frustrumRadius1^2 + ...
                    frustrumRadius1*frustrumRadius2 + ...
                    frustrumRadius2^2);
                % ...washer plus...
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(24/58*obj.total_travel_length))^2) * ...
                    (35/58-24/58) * obj.total_travel_length;
                % ...conic frustrum plus...
                frustrumRadius1 = obj.piston_diam/2 + ...
                    obj.gapProfile(35/58 * obj.total_travel_length);
                frustrumRadius2 = obj.piston_diam/2;
                V = V + pi/3 * (46/58-35/58) * obj.total_travel_length ...
                    * (frustrumRadius1^2 + ...
                    frustrumRadius1 * frustrumRadius2 + ...
                    frustrumRadius2^2);
                % ...washer
                V = V + pi * (...
                    (obj.piston_diam/2 + obj.gapProfile(58/58*obj.total_travel_length))^2) * ...
                    (x-46/58) * obj.total_travel_length;
                
            elseif x < 0 || x > 1
                V = 0;
                warning('Queried geometry outside of bounds')
            end
        end
        
        function plotGapLength(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.gapProfile(xVec(i));
            end
            plot(xVec, yVec);
        end
        
        function plotVolumeFn(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.rearVolume(xVec(i));
            end
            plot(xVec, yVec, '.');
        end
        
        function plotGapArea(obj)
            xVec = linspace(0, obj.total_travel_length, 10000);
            yVec = nan(size(xVec));
            for i = 1:length(xVec)
                yVec(i) = obj.gapArea(xVec(i));
            end
            plot(xVec, yVec, '.');
        end
    end
end