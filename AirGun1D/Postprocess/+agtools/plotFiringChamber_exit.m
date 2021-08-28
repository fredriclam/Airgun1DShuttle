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

function plotFiringChamber_exit(fullState)

% Pull model data
t = [fullState.t];
pS = [fullState.portStates];
massRate = [pS.massFlowPort];

plot(1e3*t, massRate, 'k', 'LineWidth', 1);

xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');