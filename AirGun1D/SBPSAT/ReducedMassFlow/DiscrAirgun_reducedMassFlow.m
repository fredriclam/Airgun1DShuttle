classdef DiscrAirgun_reducedMassFlow < noname.Discretization
    properties
        name        = 'Euler Wall conditions' %Short description
        description = '1D Euler equation with wall boundary conditions on both sides.' %Longer description
        order        %Order of accuracy
        schm

        physConst

        RHS

        t0
        q0
        bubble0

        H
        h
        bc_data_l
        bc_data_r
    end

    methods
        function obj = DiscrAirgun_reducedMassFlow(m,order,airgunPressure,airgunLength,airgunPortArea,airgunDepth, Mfac)
            default_arg('m',100)
            default_arg('order',4)

            %[physConst, t0, icAirgun, icBubble] = configAirgun('shocktube');
            %[physConst, t0, icAirgun, icBubble] = configAirgun('nonphys');
            [physConst, t0, icAirgun, icBubble] = configAirgun('Bolt1500LL',airgunPressure,airgunLength,airgunPortArea,airgunDepth);
            
            A  = physConst.A;
            xlim = {-physConst.L, 0};

            fprintf('Starting at time t0 = %f\n',t0)


            m = (xlim{2}-xlim{1})*m+1;

            schm = scheme.Euler1d(m,xlim,order,[],[],true);
  
            closure_l = schm.boundary_condition('l', 'wall');
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');  %% should set pressure
            L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];

            function [dq, dBubble] = RHS(q,t,bubble)
                flowState = schm.flowStateR(q);

                if t >= physConst.AirgunCutoffTime || flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    dq = q.*0;
                    dBubble = bubbleRHS_reducedMassFlow(bubble, 0, 0, 0, 0, 0, physConst, Mfac);
                    return
                end

                p_b = bubblePressure(bubble, physConst);
                dq = schm.D(q) + closure_l(q);


                %% Apply conditions on airgun and calculate q_hat
                if flowState == scheme.Euler1d.SUBSONIC_OUTFLOW

                    dq = dq + closure_r_out_sub(q,p_b);

                    % Set
                    pIn = [3];
                    pOut = [1 2];
                    permutation = [pIn pOut];
                    invPermutation(permutation) = 1:3;

                    % q_r = schm.e_R'*q;
                    % w = inv(schm.T(q_r))*q_r;

                    % Ltemp = L(q_r(1), q_r(2)/q_r(1), q_r(3));
                    % g = p_b;
                    % T = schm.T(q_r);

                    % wHat = zeros(3,1);
                    % wHat(pIn) = inv(Ltemp*T(:,pIn))*(g - Ltemp*T(:,pOut)*w(pOut)); % g := p_b
                    % wHat(pOut) = w(pOut);

                    % qHat = T*wHat;
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    % No bc required
                    qHat = schm.e_R'*q;
                else
                    % Will not run due to earlier if
                    q_r = schm.e_R'*q;
                    printExpr('q_r(1)');
                    printExpr('q_r(2)');
                    printExpr('q_r(3)');

                    % q
                    c = schm.c(q_r)
                    v = q_r(2)/q_r(1)
                    error('Undefined behaviour. flowState = %d, t = %f', flowState, t)
                end

                qHat = schm.e_R'*q;
                rho_a = qHat(1);
                v_a = qHat(2)/qHat(1);
                e_a = qHat(3);
                p_a = schm.p(qHat);
                dBubble = bubbleRHS_reducedMassFlow(bubble, rho_a, v_a, e_a, p_a, A, physConst, Mfac);
                %dBubble = bubbleRHS(bubble, rho_a, v_a, e_a, p_a, A, physConst);
                
%                 % turbulent energy dissipation
%                 turbCorrecCoeff = C; %coefficient of turbulent energy dissipation
%                 gun_area = 12.5;
%                 port_area = 12;
%                 deltaP = C*rho_a*v_a^2*(gun_area/port_area)^2;
%                 %dBubble = bubbleRHS_turbCorrec(bubble, rho_a, v_a, e_a, p_a - deltaP, A, physConst);
%                 dBubble = bubbleRHS_turbCorrec(bubble, rho_a, v_a, e_a, p_a, A, physConst, dissp_boolean, beta);
            end

            H = schm.H;

            obj.t0 = t0;

            rho0 = icAirgun.rho0(schm.u);
            rv0 = icAirgun.rv0(schm.u);
            e0 = icAirgun.e0(schm.u);

            obj.q0 = reshape([rho0'; rv0'; e0'],3*m,1);


            obj.bubble0 = [
                icBubble.R;
                icBubble.Rdot;
                icBubble.m;
                icBubble.E;
            ];

            obj.order = order;
            obj.RHS = @RHS;
            obj.H = H;
            obj.schm = schm;
            obj.h = schm.h;
            obj.physConst = physConst;
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

        function opt = defaultOpt(~, opt)
            default_arg('opt',[]);

            default_field(opt, 'method', 'rk4');
            default_field(opt, 'k', []);
            default_field(opt, 'cfl', []);
        end

        function ts = getTimestepper(obj, opt)
            default_arg('opt',[]);
            opt = obj.defaultOpt(opt);

            if isempty(opt.k)
                k = obj.getTimestep(opt);
            else
                k = opt.k;
            end

            t0 = obj.t0;
            v0 = [obj.q0; obj.bubble0];

            Nq = length(obj.q0);

            function dv = F(v,t)
                q = v(1:Nq);
                bubble = v(Nq+1:end);

                [dq, dBubble] = obj.RHS(q,t,bubble);
                dv = [dq; dBubble];
            end


            switch opt.method
                case 'rk4'
                    ts = time.Rungekutta4proper(@F,k,t0,v0);
                case 'ode15s'
                    error()
                % case 'rk4euler'
                %     t = 0;
                %     % In this case k is not the true timestep.
                %     % It will be scaled by the largest e.value.
                %     ts = Rk4euler(obj.RHS,k,t,obj.q0,@(Q)obj.schm.c(Q));
                otherwise
                    error('No timestepper method ''%s''',method)
            end
        end

        function k = getTimestep(obj, opt)
            default_arg('opt',[]);
            opt = obj.defaultOpt(opt);

            % Define table of default cfls
            cflConst = Dictionary();
            cflConst('rk4',3) = 1e-3;
            cflConst('rk4',4) = 6.0026e-01;
            cflConst('rk4',6) = 7.7568e-01;
            cflConst('rk4',8) = 9.2384e-01;

            cflConst('rk4euler',3) = 0.61829/10;
            cflConst('rk4euler',4) = 0.33078/20;
            cflConst('rk4euler',6) = 0.42126/20;
            cflConst('rk4euler',8) = 0.49837/20;

            % If no cfl was provided use a table value
            if isempty(opt.cfl)
                cfl = cflConst(opt.method, obj.order);
                cfl = cfl/2; %% For safety
            end

            k = cfl * obj.schm.h;
        end

        function repr = getTimeSnapshot(obj,ts)
            default_arg('ts',0);

            repr.airgun.x = obj.schm.x;
            repr.airgun.h = obj.schm.h;

            if ts == 0
                repr.t = obj.t0;
                repr.q  = obj.q0;
                repr.bubble = obj.bubble0;
                repr.bc = [];
                return
            end

            %% We now expect ts to be a proper timestepper

            repr.t = ts.t;

            % Extract bubble and airgun part
            v = ts.getV();
            repr.q = v(1:length(obj.q0));
            repr.bubble = v(length(obj.q0)+1:end);

            % Determine what BC is active
            repr.bc.r = cell(1,4);
            flowState = obj.schm.flowStateR(repr.q);
            if flowState == scheme.Euler1d.SUBSONIC_OUTFLOW
                repr.bc.r{4} = bubblePressure(repr.bubble, obj.physConst);
            elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                % No conditions.
            end
        end

        function [update, fh] = setupPlot(obj,~)
            lim.rho = [0 150];
            v_max = 400;
            lim.v = [-v_max v_max];
            lim.e = [0 3.2e7];
            lim.p = [0 1.2e7];

            lim.rho = [];
            lim.v = [];
            lim.e = [];
            lim.p = [];

            fh = figure();
            fh.Units = 'inches';
            fh.Position(3:4) = [16 8];
            [updateAirgun, axAirgun] = setupEulerPlot(obj,lim);

            for i = 1:length(axAirgun)
                axAirgun(i).Position(1) = -10;  % subplots delete anything they overlap so move for now
                axAirgun(i).Position(3) = 0.334659;
            end


            [updateBubble, axBubble] = setupBubbleTimeseries(obj.physConst);
            for i = 1:length(axBubble)
                axBubble(i).Position(1) = 0.570341;
                axBubble(i).Position(3) = 0.334659;
            end

            for i = 1:length(axAirgun)
                axAirgun(i).Position(1) = 0.130000; %% Restore position.
            end

            p = @(q)obj.schm.p(q);
            fhQ = figure();
            fhQ.Units = 'inches';
            fhQ.Position(3:4) = [8 8];
            updateQ = setupQTimeseries(p,'Q',@(q)obj.schm.c(q));

            fh.Position(1) = fh.Position(1) - 8;
            fhQ.Position(1) = fh.Position(1) + 16;

            e_R = obj.schm.e_R;
            function update_fun(repr)
                t = repr.t;
                q = repr.q;
                updateAirgun(t, q, repr.bc);
                updateBubble(t, repr.bubble);
                updateQ(t,e_R'*q);
            end
            repr = obj.getTimeSnapshot(0);
            update_fun(repr);

            update = @update_fun;
        end

    end

    methods(Static)
    end
end

function [update] = setupQTimeseries(p, titleStr,c)
    subplot(4,1,1)
    title(titleStr)
    [update_rho, lh_rho] = anim.setup_time_quantity_plot();
    % lh_R.Marker = '.';
    lh_rho.LineWidth = 2;
    lh_rho.Color = Color.blue;

    ylabel('rho')
    % ylim([0 10])

    subplot(4,1,2)
    [update_v, lh_v] = anim.setup_time_quantity_plot();
    [update_c, lh_c] = anim.setup_time_quantity_plot();
    % lh_m.Marker = '.';
    lh_v.LineWidth = 2;
    lh_v.Color = Color.blue;
    lh_c.LineWidth = 1;
    lh_c.Color = Color.yellow;
    ylabel('v')
    % ylim([0 10])

    subplot(4,1,3)
    [update_e, lh_e] = anim.setup_time_quantity_plot();
    % lh_E.Marker = '.';
    lh_e.LineWidth = 2;
    lh_e.Color = Color.blue;
    ylabel('e')
    % ylim([0 10])

    subplot(4,1,4)
    [update_p, lh_p] = anim.setup_time_quantity_plot();
    % lh_E.Marker = '.';
    lh_p.LineWidth = 2;
    lh_p.Color = Color.blue;
    ylabel('p')
    % ylim([0 10])

    xlabel('t')

    function update_fun(t, q)
        update_rho(t, q(1));
        update_v(t, q(2)/q(1));
        update_c(t, c(q));
        update_e(t, q(3));
        update_p(t, p(q));
    end

    update = @update_fun;
end

function [update, ax] = setupBubbleTimeseries(physConst)

    ax(1) = subplot(4,1,1);
    title('Bubble')
    [update_R, lh_R] = anim.setup_time_quantity_plot();
    % lh_R.Marker = '.';
    lh_R.LineWidth = 2;
    lh_R.Color = Color.blue;

    ylabel('R')
    % ylim([0 10])

    ax(2) = subplot(4,1,2);
    [update_p, lh_p] = anim.setup_time_quantity_plot();
    % lh_p.Marker = '.';
    lh_p.LineWidth = 2;
    lh_p.Color = Color.blue;

    ylabel('p')
    % ylim([0 10])

    ax(3) = subplot(4,1,3);
    [update_m, lh_m] = anim.setup_time_quantity_plot();
    % lh_m.Marker = '.';
    lh_m.LineWidth = 2;
    lh_m.Color = Color.blue;
    ylabel('m')
    % ylim([0 10])

    ax(4) = subplot(4,1,4);
    [update_E, lh_E] = anim.setup_time_quantity_plot();
    % lh_E.Marker = '.';
    lh_E.LineWidth = 2;
    lh_E.Color = Color.blue;
    ylabel('E')
    % ylim([0 10])

    xlabel('t')

    function update_fun(t, bubble)
        update_R(t, bubble(1));
        update_p(t, bubblePressure(bubble,physConst));
        update_m(t, bubble(3));
        update_E(t, bubble(4));

    end

    update = @update_fun;
end
