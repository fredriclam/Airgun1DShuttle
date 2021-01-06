function [update, axis_handles] = setupEulerPlot(discr, lim)
    x = discr.schm.u;

    axis_handles(1) = subplot(4,1,1);
    [update_rho, plot_handles] = anim.setup_1d_plot(x);
    plot_handles(1).LineWidth = 2;
    plot_handles(1).Color = Color.blue;
    % showGrid(lim.rho(1));
    showGrid(0);
    ylabel('rho');


    axis_handles(2) = subplot(4,1,2);
    [update_v, plot_handles] = anim.setup_1d_plot(x,{@(v,~)v, @(~,c)c, @(~,c)-c},false);
    plot_handles(1).LineWidth = 2;
    plot_handles(1).Color = Color.blue;
    plot_handles(2).LineWidth = 1;
    plot_handles(2).Color = Color.yellow;
    plot_handles(3).LineWidth = 1;
    plot_handles(3).Color = Color.yellow;
    % showGrid(lim.v(1));
    showGrid(0);
    ylabel('v');


    axis_handles(3) = subplot(4,1,3);
    [update_e, plot_handles] = anim.setup_1d_plot(x,[],false);
    plot_handles(1).LineWidth = 2;
    plot_handles(1).Color = Color.blue;
    % showGrid(lim.e(1));
    showGrid(0);
    ylabel('e');


    axis_handles(4) = subplot(4,1,4);
    [update_p, plot_handles] = anim.setup_1d_plot(x,[],false);
    plot_handles(1).LineWidth = 2;
    plot_handles(1).Color = Color.blue;
    ylabel('p');
    % showGrid(lim.p(1));
    showGrid(0);

    xlabel('x');




    update_flowstate_l = showFlowstate(axis_handles,'l');
    update_flowstate_r = showFlowstate(axis_handles,'r');

    update_bc_data_l = showBC(axis_handles,'l');
    update_bc_data_r = showBC(axis_handles,'r');

    % uistack(plot_handles(1),'top')

    did_warn = false;
    function update_fun(t,q,bc)
        default_arg('bc', struct());

        default_field(bc,'l', cell(1,4));
        default_field(bc,'r', cell(1,4));

        Q = reshape(q,3,discr.schm.m);

        if ~isreal(q) && ~did_warn
            fprintf('\nThe solution is not real!\n\n');
            did_warn = true;
        end

        rho = Q(1,:);
        v   = Q(2,:)./Q(1,:);
        e   = Q(3,:);
        p   = discr.schm.p(Q);
        c   = discr.schm.c(Q);

        update_rho(t,real(rho));
        update_v(t,real(v),real(c));
        update_e(t,real(e));
        update_p(t,real(p));

        % c = [discr.schm.c(Q(:,[1 end])), discr.schm.c(Q(:,end))];
        update_flowstate_l(v(1)  , discr.schm.flowStateL(q));
        update_flowstate_r(v(end), discr.schm.flowStateR(q));

        update_bc_data_l(bc.l);
        update_bc_data_r(bc.r);
    end

    update = @update_fun;

    function showGrid(y)
        grid_handle = line(x,x*0 + y);
        grid_handle.LineStyle = 'none';
        grid_handle.Marker = '.';
        grid_handle.Color = Color.red;
        grid_handle.MarkerSize = 6;
    end
end
function update = showFlowstate(axes_hand, boundary)
    % Constants for positioning the markers within the axes.
    % Relative coords with origin in lower left corner.
    pos_x = 0.05;
    pos_y = 0.8;

    switch boundary
        case 'l'
            % do nothing
        case 'r'
            pos_x = 1 - pos_x;
        otherwise
            error(boundary)
    end

    % Setup the markers
    for i = 1:length(axes_hand)
        axes(axes_hand(i));

        % Calculate position
        xlim = axes_hand(i).XLim;
        ylim = axes_hand(i).YLim;
        x = xlim(1) + pos_x*(xlim(2) - xlim(1));
        y = ylim(1) + pos_y*(ylim(2) - ylim(1));

        % Add the line object
        ha(i) = line(x, y);
        ha(i).LineStyle = 'none';
        ha(i).Marker = 'o';
        ha(i).MarkerSize = 12;
        ha(i).MarkerFaceColor = Color.green;
        ha(i).MarkerEdgeColor = 'none';
    end

    function update_fun(v, flowstate);
        eps = 10e-5;
        for i = 1:length(ha)
            if v > eps
                ha(i).Marker = '>';
            elseif v < -eps
                ha(i).Marker = '<';
            else
                ha(i).Marker = 'o';
            end

            if flowstate == scheme.Euler1d.SUPERSONIC_INFLOW || flowstate == scheme.Euler1d.SUPERSONIC_OUTFLOW
                ha(i).MarkerFaceColor = Color.red;
            else
                ha(i).MarkerFaceColor = Color.green;
            end
        end
    end
    update = @update_fun;
end

% g is a 1x4 cell array containging either
% an empty vector or bc data function for {rho, v, e, p};
function update = showBC(axesHand, boundary)
    xlim = axesHand.XLim;
    switch boundary
        case 'l'
            x = xlim(1);
        case 'r'
            x = xlim(2);
        otherwise
            error(boundary)
    end

    lh = cell(1,4);
    function update_fun(bc)
        for i = 1:length(bc)
            if ishandle(lh{i})
                delete(lh{i});
            end

            if isempty(bc{i})
                continue
            end

            axes(axesHand(i));
            lh{i} = line(x,bc{i});
            lh{i}.LineStyle = 'none';
            lh{i}.Marker = '.';
            lh{i}.MarkerSize = 24;
            lh{i}.Color = Color.red;
        end
    end
    update = @update_fun;
end