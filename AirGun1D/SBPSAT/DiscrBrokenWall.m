classdef DiscrBrokenWall < noname.Discretization
    properties
        name        = 'Euler Wall conditions' %Short description
        description = '1D Euler equation with wall boundary conditions on both sides.' %Longer description
        order        %Order of accuracy
        schm

        RHS

        H
        q0
        h

        bc_l
        bc_r_in
        bc_r_out
    end

    methods
        function obj = DiscrBrokenWall(m, order)
            default_arg('m',100)
            default_arg('order',4)


            p_inside = 0.35;
            p_outside = 0.2;
            rho_outside = 0.4;
            e_outside = 0.7;

            ic = icConstant(p_inside);

            xlim = {-1, 1};

            m = (xlim{2}-xlim{1})*m+1;

            rho0_fun = ic.rho0;
            rv0_fun  = ic.rv0;
            e0_fun   = ic.e0;

            schm = scheme.Euler1d(m,xlim,order,[],[],true);

            c = 1;
            gamma = 1.4;
            x_l = xlim{1};
            x_r = xlim{2};

            q_l = [rho0_fun(x_l); rv0_fun(x_l); e0_fun(x_l)];
            q_r = [rho0_fun(x_r); rv0_fun(x_r); e0_fun(x_r)];

            rho_l = @(t)q_l(1);
            v_l = @(t)q_l(2)/q_l(1);
            p_l = @(t)schm.p(q_l);
            rho_r = @(t)q_r(1);
            v_r = @(t)q_r(2)/q_r(1);
            p_r = @(t)schm.p(q_r);


            Ti_l = inv(schm.T(q_l));
            Ti_r = inv(schm.T(q_r));
            w_data_l = @(t)Ti_l*q_l;
            w_data_r = @(t)Ti_r*q_r;


            closure_l = schm.boundary_condition('l','wall');
            closure_r_wall = schm.boundary_condition('r','wall');
            closure_r_out = schm.boundary_condition('r','outflow');

            L_r_in = @(rho,u,e)[
                1 0 0;
                0 0 1;
            ];
            g_r_in = [rho_outside; e_outside];
            closure_r_in = schm.boundary_condition_L('r', L_r_in, [1 3]);

            obj.bc_l = {[], 0, [], []};
            obj.bc_r_in  = {rho_outside, [], e_outside, []};
            obj.bc_r_out = {[], [], [], p_outside};

            function dq = RHS(q,t)
                % dq = schm.D(q) + closure_l(q,t) + closure_r(q,t);

                dq = schm.D(q) + closure_l(q);


                if t < 0.15
                    dq = dq + closure_r_wall(q);
                else
                    flowState = schm.flowStateR(q);
                    if flowState == scheme.Euler1d.SUPERSONIC_INFLOW
                        warning('No bc for supersonic inflow')
                    elseif flowState == scheme.Euler1d.SUBSONIC_INFLOW
                        dq = dq + closure_r_in(q, g_r_in);
                    elseif flowState == scheme.Euler1d.SUBSONIC_OUTFLOW
                        dq = dq + closure_r_out(q, p_outside);
                    else
                        % No bc allowed for supersonic outflow.
                    end
                end
            end

            H = schm.H;

            rho0 = rho0_fun(schm.u);
            rv0 = rv0_fun(schm.u);
            e0 = e0_fun(schm.u);

            obj.q0 = reshape([rho0'; rv0'; e0'],3*m,1);

            obj.order = order;
            obj.RHS = @RHS;
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
            default_arg('method', 'rk4euler');
            default_arg('opt',struct());
            default_field(opt, 'k',[]);
            default_field(opt, 'cfl', []);

            if isempty(opt.k)
                k = obj.getTimestep(method,opt.cfl);
            end

            switch method
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

            cflConst('rk4euler',3) = 9.8608e-01;
            cflConst('rk4euler',4) = 5.0002e-01;
            cflConst('rk4euler',6) = 5.3354e-01;
            cflConst('rk4euler',8) = 2.6597e-01;

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
                repr.bc = [];
                return
            end

            % We now expect ts to be a proper timestepper
            repr.t = ts.t;
            repr.q = ts.getV();

            repr.bc.l = cell(1,4);
            repr.bc.l{2} = 0;

            if ts.t < 0.15
                repr.bc.r = cell(1,4);
                repr.bc.r{2} = 0;
            else
                flowState = obj.schm.flowStateR(repr.q);
                if flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    repr.bc.r = cell(1,4);
                elseif flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    repr.bc.r = obj.bc_r_in;
                else %Subsonic outflow
                    repr.bc.r = obj.bc_r_out;
                end
            end
        end

        function saveFrame = setupMov(obj, file)
            % if makemovies
                % save_frame = anim.setup_fig_mov(figure_handle,dirname);
            % end
            saveFrame = @()(10);
        end

        function [update, figure_handle] = setupPlot(obj,~)
            lim.rho = [-0.1 1.5];
            lim.v   = [-1.1 2.3];
            lim.e   = [-0.1 2.1];
            lim.p   = [-0.1 1.1];

            [update_plots, figure_handle] = setupEulerPlot(obj,lim);

            function update_fun(repr)
                update_plots(repr.t, repr.q, repr.bc);
            end
            repr = obj.getTimeSnapshot(0);
            update_fun(repr)

            update = @update_fun;
        end

    end

    methods(Static)
    end
end