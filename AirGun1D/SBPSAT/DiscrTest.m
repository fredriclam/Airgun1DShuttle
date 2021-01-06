classdef DiscrTest < noname.Discretization
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
        function obj = DiscrTest(m, order, ic)
            default_arg('m',100)
            default_arg('order',4)

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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% With background flow to the right
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % closure_l = schm.boundary_condition('l','char',w_data_l);
            % closure_r = schm.boundary_condition('r','outflow_rho',rho_r);
            % obj.bc_data_l = {};
            % obj.bc_data_r = {rho_r, [], [], []};

            % closure_l = schm.boundary_condition('l','inflow_rho',rho_l,v_l);
            % closure_r = schm.boundary_condition('r','char',w_data_r);
            % obj.bc_data_l = {rho_l, v_l, [],[]};
            % obj.bc_data_r = {};

            % closure_l = schm.boundary_condition('l','inflow_rho',rho_l,v_l);
            % closure_r = schm.boundary_condition('r','outflow_rho',rho_r);
            % obj.bc_data_l = {rho_l, v_l, [],[]};
            % obj.bc_data_r = {rho_r};

            % closure_l = schm.boundary_condition('l','inflow',p_l,v_l);
            % closure_r = schm.boundary_condition('r','outflow',p_r);
            % obj.bc_data_l = {[], v_l, [], p_l};
            % obj.bc_data_r = {[], [], [], p_r};


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% With background flow to the left
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % closure_l = schm.boundary_condition('l','char', w_data_l);
            % closure_r = schm.boundary_condition('r','inflow_rho',rho_r,v_r);
            % obj.bc_data_l = {};
            % obj.bc_data_r = {rho_r, v_r, [],[]};

            % closure_l = schm.boundary_condition('l','outflow_rho', rho_l);
            % closure_r = schm.boundary_condition('r','char', w_data_r);
            % obj.bc_data_l = {rho_l, [], [], []};
            % obj.bc_data_r = {};

            % closure_l = schm.boundary_condition('l','outflow_rho',rho_l);
            % closure_r = schm.boundary_condition('r','inflow_rho',rho_r,v_r);
            % obj.bc_data_l = {rho_l};
            % obj.bc_data_r = {rho_r, v_r, [],[]};






            % closure_l_in  = schm.boundary_condition('l','inflow' ,p_l,v_l);
            % closure_l_out = schm.boundary_condition('l','outflow',p_l);
            % closure_r_in  = schm.boundary_condition('r','inflow' ,p_r,v_r);
            % closure_r_out = schm.boundary_condition('r','outflow',p_r);

            closure_l_in  = schm.boundary_condition('l','inflow_rho' ,rho_l,v_l);
            closure_l_out = schm.boundary_condition('l','outflow_rho',rho_l);
            closure_r_in  = schm.boundary_condition('r','inflow_rho' ,rho_r,v_r);
            closure_r_out = schm.boundary_condition('r','outflow_rho',rho_r);

            closure_l_char = schm.boundary_condition('l','char', w_data_l);
            closure_r_char = schm.boundary_condition('r','char', w_data_r);
            inwin = @(t)step(t, 0.3, 0.35);
            function dq = RHS(q,t)
                % dq = schm.D(q) + closure_l(q,t) + closure_r(q,t);

                dq = schm.D(q);

                switch schm.flowStateL(q)
                    case scheme.Euler1d.SUBSONIC_INFLOW
                        dq = dq + closure_l_in(q,t);
                    case scheme.Euler1d.SUBSONIC_OUTFLOW
                        dq = dq + closure_l_out(q,t);
                    otherwise
                        error('Flowstate not handled: %d', schm.flowStateL(q));
                end

                switch schm.flowStateR(q)
                    case scheme.Euler1d.SUBSONIC_INFLOW
                        dq = dq + closure_r_in(q,t);
                    case scheme.Euler1d.SUBSONIC_OUTFLOW
                        dq = dq + closure_r_out(q,t);
                    otherwise
                        error('Flowstate not handled: %d', schm.flowStateR(q));
                end

                % if t >= 1 && t < 0.35
                %     % dq = dq + closure_l_infl(q,t);
                %     dq = dq + closure_l_otfl(q,t);
                % else
                %     dq = dq + closure_l_wall(q);
                % end
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
            default_arg('k',[]);
            default_arg('cfl', []);

            if isempty(k)
                k = obj.getTimestep(method,cfl);
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

            cflConst('rk4euler',3) = 0.01829;
            cflConst('rk4euler',4) = 0.33078*0.5*0.5;
            cflConst('rk4euler',6) = 0.42126;
            cflConst('rk4euler',8) = 0.49837;

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
            lim.v   = [-1.1 2.3];
            lim.e   = [-0.1 2.1];
            lim.p   = [-0.1 1.1];

            [[update_plots, figure_handle] = setupEulerPlot(obj,lim);

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