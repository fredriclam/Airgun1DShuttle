%% Plot firing chamber wall state as a function of t
% Plots the specified state at the closed end
% Inputs:
%   state: Full state struct-vector (from compute_full_state)
%   varname: name of field (p, u, T, ...) to plot
%   (optional) rawData: struct(t:time axis, y:data axis) vector with
%   field data objects
% Outputs:
%   t: Extracted t-axis for model
%   field_L: Extracted field at wall for model
% Additional:
%   Plots to current figure.

function [C, h] = plotFiringChamber_xt_nofill( ...
    fullState, metadata, caseKeyContext, linecolor)

% Pull model data
t = [fullState.t];
eDS = [fullState.eulerDomainStates];
p = [eDS.p];
x = metadata.discretization.schm.u;

has_shuttle_position = false;

if has_shuttle_position
    tL = tiledlayout(4,4);
    nexttile(tL, 5, [3,3]);
end

[tGrid, xGrid] = meshgrid(t, x);
[C, h] = contour(xGrid(:,1:15000), 1e3*tGrid(:,1:15000), 1e-6 * p(:,1:15000), ...
    'LevelStep', 0.5, 'LineColor', linecolor, 'LineWidth', 1);

xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
% cbh = colorbar('TickLabelInterpreter', 'latex', 'FontSize', 14, ...
%     'location', 'northoutside');
% xlabel(cbh, '$p$ [MPa]', 'FontSize', 14, 'Interpreter', 'latex')
% colormap cool
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on' ...
    ...,'YTickLabel', {'','','','',''}
    );

% % Plot tracker line for next figure
% tIndexSample = 5000;
% colorTrackerLine = [252, 111, 23]/255;
% hold on
% plot([x(1) x(end)], t(5000)*[1 1], ...
%     'Color', colorTrackerLine, 'LineWidth', 0.5)
% hold off 
% ylimCurrent = ylim;
% caxis([0,7])

  
% Dummy lines for legend
% for i = 1:6
%     plot(shuttle_position(1), t(1), ...
%         [caseKeyContext.colorMap{i}, '-'])
%     hold on
% end

% Shuttle position attachment
if has_shuttle_position
    nexttile(tL, 8, [3,1]);
    for i = 1:length(caseKeyContext.caseKeySwitchIndices)-1
        % Include the right boundary element too for continuity
        plotRange = caseKeyContext.caseKeySwitchIndices(i)+1: ...
            min(length(t), caseKeyContext.caseKeySwitchIndices(i+1)+1);
    %     plot(shuttle_position(plotRange), t(plotRange), ...
    %         [caseKeyContext.colorMap{1+caseKeyContext.caseKeyHistory(plotRange(1))}], 'LineWidth', 1);
        plot(shuttle_position(plotRange), 1e3*t(plotRange), ...
            'k', 'LineWidth', 1);
        hold on
    end
    hold off
end

% legendLabels = {'Closed', ...
%     'Subsonic', ...
%     'Port choked', ...
%     'Chamber choked', ...
%     'Chamber choked*', ...
%     'Relaxation'};
% legend(legendLabels, 'Interpreter', 'latex', 'location', 'eastoutside', ...
%     'FontSize', 10)

% set(gca, 'YTickLabel', {'','','','',''}, ...
%     'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
%     'XMinorTick', 'on', 'YMinorTick', 'on')
% Sync vertical axis
% ylim([0, ylimCurrent(2)])
% set(gca, 'position', [0.76, 0.13, 0.15, 0.8])
% xlabel('$\xi$ [m]', 'Interpreter', 'latex', 'FontSize', 14)

% Resize axes
% subplot(1,3,3);
% set(gca, 'position', [0.68, 0.13, 0.15, 0.68])
% subplot(1,3,1);
% ylim([0, ylimCurrent(2)])
% set(gca, 'position', [0.07, 0.13, 0.16, 0.68])
% subplot(1,3,2);
% set(gca, 'position', [0.25, 0.13, 0.40, 0.68])
% set(gcf, 'position', [200 399 850 483])