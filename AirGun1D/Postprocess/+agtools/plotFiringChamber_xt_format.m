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

function plotFiringChamber_xt_format( ...
    fullState, metadata)

% Pull model data
t = [fullState.t];
eDS = [fullState.eulerDomainStates];
p = [eDS.p];
x = metadata.discretization.schm.u;

sS = [fullState.shuttleStates];
if isfield(sS, 'shuttle_position')
    shuttle_position = [sS.shuttle_position];
    has_shuttle_position = true;
else
    has_shuttle_position = false;
end

if has_shuttle_position
    tL = tiledlayout(4,5);
    nexttile(tL, 6, [3,3]);
end

[tGrid, xGrid] = meshgrid(t, x);
contourf(xGrid(:,1:15000), 1e3*tGrid(:,1:15000), 1e-6 * p(:,1:15000), ...
    'LevelStep', 0.15);

xlabel('{\it{x}} (m)', 'FontSize', 14)
ylabel('{\it{t}} (ms)', 'FontSize', 14)
cbh = colorbar('FontSize', 14, ...
    'location', 'northoutside');
xlabel(cbh, '{\it{p}} (MPa)', 'FontSize', 14)
colormap cool
set(gca, 'FontSize', 14, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', [0.01, 0.025]*1.6...,'YTickLabel', {'','','','',''}
    );
ylim([0, 100])

% Plot tracker line for next figure
tIndexSample = 5000;
colorTrackerLine = [252, 111, 23]/255;
hold on
plot([x(1) x(end)], t(5000)*[1 1], ...
    'Color', colorTrackerLine, 'LineWidth', 0.5)
hold off 
ylimCurrent = ylim;
caxis([0,7])

  
% Dummy lines for legend
% for i = 1:6
%     plot(shuttle_position(1), t(1), ...
%         [caseKeyContext.colorMap{i}, '-'])
%     hold on
% end

if has_shuttle_position
    % Plot shuttle position
    nexttile(tL, 9, [3,1]);
    plot(shuttle_position, 1e3*t, ...
        'k', 'LineWidth', 1);
    set(gca, 'YTickLabel', {'','','','',''}, ...
        'FontSize', 14, ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.01, 0.025]*1.6)
    % Sync vertical axis
    xlim([0, 0.07])
    ylim([0, ylimCurrent(2)])
    % set(gca, 'position', [0.76, 0.13, 0.15, 0.8])
    xlabel('\xi (m)', 'FontSize', 14)
    
    % Mass rate out
    nexttile(tL, 10, [3,1]);
    t = [fullState.t];
    pS = [fullState.portStates];
    massRate = [pS.massFlowPort];

    plot(massRate, 1e3*t, 'k', 'LineWidth', 1);

    xlabel('$\dot{\it{m}}$ (kg/s)', 'FontSize', 14, 'Interpreter', 'latex')
%     ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca, 'YAxisLocation','right', ...
        ...'YTickLabel', {'','','','',''}, ...
        'FontSize', 14, ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.01, 0.025]*1.6)
    % Sync vertical axis
    ylim([0, ylimCurrent(2)])
    ylabel("{\it{t}} (ms)")
end

% legendLabels = {'Closed', ...
%     'Subsonic', ...
%     'Port choked', ...
%     'Chamber choked', ...
%     'Chamber choked*', ...
%     'Relaxation'};
% legend(legendLabels, 'Interpreter', 'latex', 'location', 'eastoutside', ...
%     'FontSize', 10)

% Resize axes
% subplot(1,3,3);
% set(gca, 'position', [0.68, 0.13, 0.15, 0.68])
% subplot(1,3,1);
% ylim([0, ylimCurrent(2)])
% set(gca, 'position', [0.07, 0.13, 0.16, 0.68])
% subplot(1,3,2);
% set(gca, 'position', [0.25, 0.13, 0.40, 0.68])
set(gcf, 'position', [200 399 850 483])