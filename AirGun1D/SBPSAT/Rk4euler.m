% Rungekutta 4 with dynamic time step.
classdef Rk4euler < time.Timestepper
    properties
        F
        k
        t
        v
        m
        n
        c
    end


    methods
        % Runs rungekutta where the timestep is scaled by the largest eigenvalue.
        function obj = Rk4euler(F, k, t0, v0,c)
            obj.F = F;
            obj.k = k;
            obj.t = t0;
            obj.v = v0;
            obj.c = c;
            obj.m = length(v0);
            obj.n = 0;
        end

        function [v,t] = getV(obj)
            v = obj.v;
            t = obj.t;
        end

        function obj = step(obj)
            q = obj.v;
            Q = reshape(q,3,obj.m/3);

            c = obj.c(Q);
            maxeig = max(abs(Q(2,:)./Q(1,:)) + c);

            k = obj.k * maxeig;

            obj.v = time.rk4.rungekutta_4(obj.v, obj.t, k, obj.F);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end


    methods (Static)
        function k = getTimeStep(lambda)
            k = rk4.get_rk4_time_step(lambda);
        end
    end

end