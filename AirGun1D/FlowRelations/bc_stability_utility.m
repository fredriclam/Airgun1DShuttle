% Polynomial mapping for Pi Hats

function [PiSolutions, functions, alt_paths] = ...
    bc_stability_utility(Pi, M, areaRatio, is_plotting)
alt_paths = struct();
functions.get_Pi_roots = @(Pi, M, areaRatio, gamma) ...
    get_Pi_roots(Pi, M, areaRatio, gamma);
functions.get_choked_pressure_range = @(areaRatio, gamma) ...
    get_choked_pressure_range(areaRatio, gamma);
if nargin == 0
    PiSolutions = [];
    return
elseif nargin < 4
    is_plotting = false;
end
% Compute all solutions Pi using polynomial solver
% Typically provides one physical and one unphysical root Pi
[PiSolutions, roots_raw, constraintFn] = ...
    get_Pi_roots(Pi, M, areaRatio);
if is_plotting
    % Create plot of constraint equation g(Pi)
    subplot(1,2,1);
    ezplot(constraintFn)
    % Create complex plane plot of candidate solutions Pi
    subplot(1,2,2);
    plot(real(roots_raw.^7), imag(roots_raw.^7), 'k.', 'MarkerSize', 15);
    xlabel 'Re(\sigma^7)'
    ylabel 'Im(\sigma^7)'
end

%% Compute the throat choked condition (subsonic to sonic inlet)
Pi_shockedNozzle_subsonic = NaN;
gamma = 1.4;
if areaRatio > 1
    if isempty(PiSolutions)
        Pi_shockedNozzle_subsonic = Pi * ( ...
            (M + 2/(gamma-1)) / (1 + 2/(gamma-1))...
        )^(2*gamma/(gamma-1));
    end
end

alt_paths.Pi_shockedNozzle_subsonic = Pi_shockedNozzle_subsonic;

end

function [Pi_max, Pi_min] = get_choked_pressure_range( ...
    areaRatio, gamma)
    % Computes the max/min pressure ratio that can be sustained in a nozzle
    % that is choked at the inlet; above max, pressure is discontinuous
    % between the exit pressure and the back pressure (reservoir pressure),
    % and below max, pressure is insufficient to choke at the inlet.
if nargin < 4
    gamma = 1.4;
end
if gamma ~= 1.4
    error("Not implemented for gamma other than 1.4 == 7/5.")
end

% Define polynomial form
poly_form = [ ...
    ((gamma-1)/2)  / areaRatio^2 ...         % 12
    0 ...                                    % 11
    0 ...                                    % 10
    0 ...                                    % 9
    0 ...                                    % 8
    0 ...                                    % 7
    0 ...                                    % 6
    0 ...                                    % 5
    0 ...                                    % 4
    0 ...                                    % 3
    -(1 + (gamma-1)/2) ...                   % 2
    0 ...                                    % 1
    1                                        % 0
];
% Compute complex polynomial roots of mapped polynomial problem
roots_raw = roots(poly_form);
Pi_min = min(roots_raw);
Pi_max = max(roots_raw);

%% Root filtering in Pi-space
% Define explicit hat constraint (admitting however many solutions)
% for evaluating candidate solutions Pi
eqn_grid2target = @(Pi) (1 + (gamma-1)/2) * Pi.^((gamma-1)/(gamma)) ...
    - (gamma-1)/2 / areaRatio.^2 .* Pi.^((gamma+1)/gamma) - 1;
PiSolutions = roots_raw.^7;
% Filter out complex solutions (imag part should be exactly zero)
PiSolutions(angle(roots_raw) ~= 0) = [];
% Filter out solutions Pi that do not satisfy the constraint explicitly
% using an arbitrary small tolerance
for i = 1:length(PiSolutions)
    if abs(eqn_grid2target(PiSolutions(i))) > 1e-5
        PiSolutions(i) = NaN;
    end
end
PiSolutions(isnan(PiSolutions)) = [];

% Get outputs
Pi_min = min(PiSolutions);
Pi_max = max(PiSolutions);

end

function [PiSolutions, roots_raw, constraintFn] = get_Pi_roots( ...
    Pi, M, areaRatio, gamma)
if nargin < 4
    gamma = 1.4;
end
if gamma ~= 1.4
    error("Not implemented for gamma other than 1.4 == 7/5.")
end
% Compute problem invariant
J = Pi.^((gamma-1)/(2*gamma)) .* (M + 2/(gamma-1));
% Define polynomial form
poly_form = [ ...
    (2 / (gamma - 1))^2 / areaRatio^2 ... % 12
    -4/(gamma-1) * J / areaRatio^2 ...    % 11
    J^2 / areaRatio^2 ...                 % 10
    0 ...                              % 9
    0 ...                              % 8
    0 ...                              % 7
    0 ...                              % 6
    0 ...                              % 5
    0 ...                              % 4
    0 ...                              % 3
    -2 * (gamma+1) / (gamma-1)^2 ...   % 2
    4 / (gamma - 1) * J ...            % 1
    -(J^2 - 2/(gamma-1))               % 0
];
% Compute complex polynomial roots of mapped polynomial problem
roots_raw = roots(poly_form);

%% Root filtering in Pi-space
% Define explicit hat constraint (admitting however many solutions)
% for evaluating candidate solutions Pi
eqn_grid2target = @(Pi) J - Pi.^((gamma-1)/(2*gamma)) ...
    .* (sqrt( ...
        2/(gamma-1)*(1-Pi.^((gamma-1)/gamma)) ...
        ./(Pi.^((gamma-1)/gamma)-Pi.^((gamma+1)/gamma) ./ areaRatio.^2)) ...
    + 2/(gamma-1));
PiSolutions = roots_raw.^7;
% Filter out complex solutions (imag part should be exactly zero)
PiSolutions(angle(roots_raw) ~= 0) = [];
% Filter out solutions Pi that do not satisfy the constraint explicitly
% using an arbitrary small tolerance
for i = 1:length(PiSolutions)
    if abs(eqn_grid2target(PiSolutions(i))) > 1e-5
        PiSolutions(i) = NaN;
    end
end
PiSolutions(isnan(PiSolutions)) = [];
% Sort output solutions
PiSolutions = sort(PiSolutions);

% Export constraint function used
constraintFn = eqn_grid2target;

end