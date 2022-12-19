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

function [t, field_L, fieldData] = plotFiringChamber_Wall( ...
    fullState, varname, rawData)

% Pull model data
t = [fullState.t];
eDS = [fullState.eulerDomainStates];
field = [eDS.(varname)];
field_L = field(1,:);

% Set scale multipliers
tScale = 1e-3;
fieldScale = 1;
if strcmpi('p', varname)
    % Set (ms, MPa) scale multipliers
    fieldScale = 1e6;
end
plot(t/tScale, field_L/fieldScale, 'k--', 'LineWidth', 1);

if nargin >= 3
    if strcmpi('p', varname)
        % Prep field data
        psiPa_conversion = 6894.75729;
        % Voltage baseline at 0 psi
        nominalminimum = 1000*psiPa_conversion; 
        % Manual data index input
        begin_index = 3700*5;
        end_index = begin_index+1600*5;
        fieldData = struct();
        fieldData.t = 1e3 * ...
            (rawData(24).iNetTimeAxisP(begin_index:end_index) - ...
             rawData(24).iNetTimeAxisP(begin_index));
        fieldData.y = 1e-6* (nominalminimum + ...
            psiPa_conversion * ...
            rawData(24).iNetCh16Data(begin_index:end_index));
    elseif strcmpi('T', varname)
        % Convert experimental data from deg F to K
        F_to_K = @(F) 5/9*(F-32) + 273.15;
        T_exper_K = F_to_K(rawData(24).iNetCh4Data);
        % Manual search index input
        begin_index = 3700;
        end_index = begin_index+600;

        fieldData = struct();
        fieldData.t = 1e3* ...
            (rawData(24).iNetTimeAxisT(begin_index:end_index) - ...
             rawData(24).iNetTimeAxisT(begin_index));
        fieldData.y = T_exper_K(begin_index:end_index);
    end
    % Plot field data
    hold on
    plot(fieldData.t, fieldData.y, 'b-', 'LineWidth', 1);
end

hold off
xlabel('$t$ [ms]', 'Interpreter', 'latex', 'FontSize', 14)
if strcmpi('p', varname)
    ylabel('$p$ [MPa]', 'Interpreter', 'latex', 'FontSize', 14)
    xlim([0, 300]);
elseif strcmpi('T', varname)
    ylabel('$T$ [K]', 'Interpreter', 'latex', 'FontSize', 14)
    xlim([0, 300]);
    ylim([0, 300]);
else
    ylabel("$" + varname + "$", 'Interpreter', 'latex', 'FontSize', 14)
end
legend({'Model', ...
        'Data', ...
        ...
       }, ...
      'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on')
end