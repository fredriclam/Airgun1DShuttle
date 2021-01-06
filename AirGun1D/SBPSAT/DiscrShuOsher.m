classdef DiscrShuOsher < noname.Discretization
    properties
        name        = 'Euler Wall conditions' %Short description
        description = '1D Euler equation with wall boundary conditions on both sides.' %Longer description
        order        %Order of accuracy
        schm

        RHS

        H
        q0
        h
        bc_data_l
        bc_data_r
    end

    methods
        function obj = DiscrShuOsher(m,order)
            default_arg('m',100)
            default_arg('order',4)

            ic = icShuOsher();

            xlim = {0, 1};

            m = (xlim{2}-xlim{1})*m+1;

            schm = scheme.Euler1d(m,xlim,order,[],[],true);

            x_l = xlim{1};
            q_l = [ic.rho0(x_l); ic.rv0(x_l); ic.e0(x_l)];
            Ti_l = inv(schm.T(q_l));
            w_data_l = @(t)Ti_l*q_l;

            rho_l = ic.rho0(x_l);
            p_l   = ic.p0(x_l);
            v_l   = ic.rv0(x_l)/ic.rho0(x_l);

            L = @(rho, u, ~)[
                1          0   0;
                0      1/rho   0;
                0 -1/2*u*0.4 0.4;
            ];

            g = [
                rho_l;
                v_l;
                p_l;
            ];

            [closure_l] = schm.boundary_condition_L('l', L, [1 2 3]);
            [closure_r] = schm.boundary_condition('r', 'wall');

            RHS = @(q,t) schm.D(q) + closure_l(q,g) + closure_r(q);

            obj.bc_data_l = {@(t)rho_l, @(t)v_l, [],@(t)p_l};
            obj.bc_data_r = {[],@(t)0,[],[]};

            H = schm.H;

            rho0 = ic.rho0(schm.u);
            rv0 = ic.rv0(schm.u);
            e0 = ic.e0(schm.u);

            obj.q0 = reshape([rho0'; rv0'; e0'],3*m,1);

            obj.order = order;
            obj.RHS = RHS;
            obj.H = H;
            obj.schm = schm;
            obj.h = schm.h;
        end

        function printInfo(obj)
            fprintf('%s\n', obj.name);
            fprintf('n = %d\n', size(obj));
            fprintf('h = %f\n', obj.schm.h);
            fprintf('xlim = (%.2f,%.2f)\n', obj.schm.x(1),obj.schm.x(end));
        end

        % Return the number of DOF
        function n = size(obj)
            n = length(obj.q0);
        end

        function ts = getTimestepper(obj, opt)
            default_arg('opt',struct());

            default_field(opt, 'method', 'rk4euler');
            default_field(opt, 'k', []);
            default_field(opt, 'cfl', []);

            if isempty(opt.k)
                k = obj.getTimestep(opt.method,opt.cfl);
            end

            switch opt.method
                case 'rk4'
                    t = 0;
                    ts = time.Rungekutta4proper(obj.RHS,k,t,obj.q0);
                case 'rk4euler'
                    t = 0;
                    % In this case k is not the true timestep.
                    % It will be scaled by the largest e.value.
                    ts = Rk4euler(obj.RHS,k,t,obj.q0,@(Q)obj.schm.c(Q));
                otherwise
                    error('No timestepper method ''%s''',method)
            end
        end

        function k = getTimestep(obj, method, cfl)
            default_arg('cfl',[])


            % Define table of default cfls
            cflConst = Dictionary();
            cflConst('rk4',3) = 1.1266;
            cflConst('rk4',4) = 6.0026e-01;
            cflConst('rk4',6) = 7.7568e-01;
            cflConst('rk4',8) = 9.2384e-01;

            cflConst('rk4euler',3) = 0.61829/10;
            cflConst('rk4euler',4) = 0.33078/20;
            cflConst('rk4euler',6) = 0.42126/20;
            cflConst('rk4euler',8) = 0.49837/20;

            % If no cfl was provided use a table value
            if isempty(cfl)
                cfl = cflConst(method,obj.order);
                cfl = cfl/2; %% For safety
            end

            k = cfl * obj.schm.h;
        end

        function repr = getTimeSnapshot(obj,ts)
            default_arg('ts',0);

            repr.x = obj.schm.x;
            repr.h = obj.schm.h;

            if ts == 0
                repr.t  = 0;
                repr.q  = obj.q0;
                return
            end

            % We now expect ts to be a proper timestepper
            repr.t = ts.t;
            repr.q = ts.getV();
        end

        function saveFrame = setupMov(obj, file)
            % if makemovies
                % save_frame = anim.setup_fig_mov(figure_handle,dirname);
            % end
            saveFrame = @()(10);
        end

        function [update, figure_handle] = setupPlot(obj,~)
            lim.rho = [-0.1 7];
            lim.v = [-4 4];
            lim.e = [-0.1 50];
            lim.p = [-0.1 15];

            [update_plots, figure_handle] = setupEulerPlot(obj,lim);

            function update_fun(repr)
                update_plots(repr.t, repr.q);
            end
            repr = obj.getTimeSnapshot(0);
            update_fun(repr)

            update = @update_fun;
        end

    end

    methods(Static)
    end
end