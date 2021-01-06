classdef DiscrWallTest < noname.Discretization
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
        function obj = DiscrWallTest(m,order,ic)
            default_arg('m',100)
            default_arg('order',4)

            xlim = {-1, 1};

            m = (xlim{2}-xlim{1})*m+1;

            rho0_fun = ic.rho0;
            rv0_fun  = ic.rv0;
            e0_fun   = ic.e0;

            schm = scheme.Euler1d(m,xlim,order,[],[],true);

            [closure_l] = schm.boundary_condition('l','wall');
            [closure_r] = schm.boundary_condition('r','wall');

            RHS = @(q,t) schm.D(q) + closure_l(q) + closure_r(q);

            obj.bc_data_l = {[],@(t)0,[],[]};
            obj.bc_data_r = {[],@(t)0,[],[]};

            H = schm.H;

            rho0 = rho0_fun(schm.u);
            rv0 = rv0_fun(schm.u);
            e0 = e0_fun(schm.u);

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

            cflConst = Dictionary();
            cflConst('rk4',3) = 1.1266;
            cflConst('rk4',4) = 6.0026e-01;
            cflConst('rk4',6) = 7.7568e-01;
            cflConst('rk4',8) = 9.2384e-01;

            cflConst('rk4euler',2) = 3.1252e-01;
            cflConst('rk4euler',4) = 4.8169e-01;
            cflConst('rk4euler',6) = 6.2622e-01;
            cflConst('rk4euler',8) = 6.9324e-01;
            cflConst('rk4euler',3) = 8.7622e-01/3;
            cflConst('rk4euler',5) = 8.7439e-01;
            cflConst('rk4euler',7) = 8.4659e-01;
            cflConst('rk4euler',9) = 7.9968e-01;

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
            lim.rho = [-0.1 1.5];
            lim.v = [-2.1 2.1];
            lim.e = [-0.1 3.1];
            lim.p = [-0.1 1.1];

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