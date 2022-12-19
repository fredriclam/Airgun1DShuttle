%% Create data
AR_vector = [0.8, 0.99, 1.01, 1.2];
N_runs = length(AR_vector);
gamma = 1.4;
clear computed_bdry;
for i = 1:length(AR_vector)
    computed_bdry(i) = compute_M_p (AR_vector(i), gamma);
end

%% Plot M
figure(51); clf
for i = 1:N_runs
    subplot(2,N_runs,i);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             computed_bdry(i).M)
    title ("subs \alpha = " + AR_vector(i));
    colorbar
    
    subplot(2,N_runs,i+N_runs);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             computed_bdry(i).M_sups)
    title ("sups \alpha = " + AR_vector(i));
    colorbar
end

%% Plot Pi
figure(52); clf
for i = 1:N_runs
    subplot(2,N_runs,i);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             computed_bdry(i).Pi)
    title ("subs \alpha = " + AR_vector(i));
    colorbar
    
    subplot(2,N_runs,i+N_runs);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             computed_bdry(i).Pi_sups)
    title ("sups \alpha = " + AR_vector(i));
    colorbar
end

%% Plot the penalty trajectory
figure(53); clf
for i = 1:N_runs
    subplot(2,N_runs,i);
    computed_bdry(i).subs_plot()
    title ("subs \alpha = " + AR_vector(i));
    
    subplot(2,N_runs,i+N_runs);
    computed_bdry(i).sups_plot()
    title ("sups \alpha = " + AR_vector(i));
end

%% Compute trajectories with fill-in values
figure(54); clf
colormap([1 1 1; 0 0 0])
sonicPressure = @(p, M) ...
    p .* ...
    ((1 + (gamma-1)/2) ./ (1 + (gamma-1)/2 * M.^2)).^(-gamma / (gamma-1));
for i = 1:N_runs
    subplot(2,N_runs,i);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             sonicPressure(computed_bdry(i).Pi, ...
             computed_bdry(i).M) > 1)
    computed_bdry(i).subs_clipped_plot()
    title ("subs \alpha = " + AR_vector(i));
    caxis([0, 1])
    
    subplot(2,N_runs,i+N_runs);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             sonicPressure(computed_bdry(i).Pi_sups, ...
                computed_bdry(i).M_sups) > 1)
    computed_bdry(i).sups_plot()
    title ("sups \alpha = " + AR_vector(i));
    caxis([0, 1])
end

%%
figure(55); clf
% colormap([1 1 1; 0 0 0])
colormap parula

sonicPressure = @(p, M) ...
    p .* ...
    ((1 + (gamma-1)/2) ./ (1 + (gamma-1)/2 * M.^2)).^(-gamma / (gamma-1));
for i = 1:N_runs
    subplot(1,N_runs,i);
    contourf(computed_bdry(i).M_grid, ...
             computed_bdry(i).p_grid, ...
             computed_bdry(i).Pi, 'LevelStep', 0.05, 'LineStyle', 'none')
%     computed_bdry(i).subs_clipped_plot()
    title ("subs \alpha = " + AR_vector(i));
    cb = colorbar;
    cb.Label.String = '\hat{Pi}';
    caxis([0, 2])
    hold on
    surf(computed_bdry(i).M_grid, ...
         computed_bdry(i).p_grid, ...
         -1 + 3*(sonicPressure(computed_bdry(i).Pi, ...
         computed_bdry(i).M) > 1), ...
         'LineStyle', 'none', 'FaceColor', 'black')
    hold off
    
    xlabel 'M_{R}'
    ylabel '\Pi_R'
end


%% Legacy
if false
    %% Create default plot
    figure(59); clf
    subplot(1,3,1);
    title ("\alpha = " + 0.8);
    exports1.create_plot();
    
    subplot(1,3,2);
    title ("\alpha = " + 1.2);
    exports2.create_plot();
    
    %%
    figure(60); clf
    subplot(2,2,1);
    contourf(exports1.Pi, 'LevelStep', 0.1)
    colorbar
    
    subplot(2,2,2);
    contourf(exports1.M, 'LevelStep', 0.1)
    colorbar
    
    subplot(2,2,4);
    M_subs = exports2.M;
    M_subs(exports2.flags ~= 1) = NaN;
    M_subs = real(M_subs);
    contourf(M_subs, 'LevelStep', 0.1)
    colorbar
    
    %%
    figure(61); clf
    contourf(exports2.M_grid, exports2.p_grid, ...
        exports2.M-exports2.M_grid);
    hold on
    contour(exports2.M_grid, exports2.p_grid, ...
        exports2.M-exports2.M_grid, [0 0], 'LineStyle', "-" ,'LineColor', [1 0 0]);
    contour(exports2.M_grid, exports2.p_grid, ...
        exports2.Pi-exports2.p_grid, [0 0], 'LineStyle', '--', 'LineColor', [0 0 0]);
    hold off
end

function exports = compute_M_p (AR, gamma, backend_solver)

if nargin < 3
    if gamma == 1.4
        % Robust solver for gamma = 1.4
        backend_solver = 'roots';
    else
        backend_solver = 'fsolve';
    end
end

hat_ranges = struct();
hat_ranges.M = linspace(0, 2, 100);
hat_ranges.p = linspace(0, 4, 200);
hat_ranges.p = hat_ranges.p(2:end);

% Create meshgrid of M, P
[meshgrids_ranges.M, meshgrids_ranges.p] = ...
    meshgrid(hat_ranges.M, hat_ranges.p);
% Compute BC invariant as function of M, P
meshgrids_ranges.data = meshgrids_ranges.p.^((gamma-1)/(2*gamma)) ...
    .* (meshgrids_ranges.M + 2/(gamma-1));

% Define nonlinear scalar constraint
eqn_grid2target = @(data, Pi, areaRatio) data - Pi.^((gamma-1)/(2*gamma)) ...
    .* (sqrt( ...
        2/(gamma-1)*(1-Pi.^((gamma-1)/gamma)) ...
        ./(Pi.^((gamma-1)/gamma)-Pi.^((gamma+1)/gamma) ./ areaRatio.^2)) ...
    + 2/(gamma-1));

% Solve nonlinear scalar constraint for P-hat
PiHats_2 = nan(size(meshgrids_ranges.data));
PiHats_sup = nan(size(meshgrids_ranges.data));
flagsPiHats_2 = nan(size(meshgrids_ranges.data));
flagsPiHats_sup = flagsPiHats_2;
for i = 1:size(meshgrids_ranges.data,1)
    for j = 1:size(meshgrids_ranges.data,2)
        if strcmpi('roots', backend_solver)
            PiCand = bc_stability_utility( ...
                meshgrids_ranges.p(i,j), ...
                meshgrids_ranges.M(i,j), ...
                AR);
            if length(PiCand) >= 1
                PiHats_sup(i,j) = PiCand(1);
                flagsPiHats_sup(i,j) = 1;
            else
                flagsPiHats_sup(i,j) = 0;
            end
            if length(PiCand) >= 2
                PiHats_2(i,j) = PiCand(2);
                flagsPiHats_2(i,j) = 1;
            else
                flagsPiHats_2(i,j) = 0;
            end
            
            if AR > 1 && isempty(PiCand)
                flagsPiHats_2(i,j) = 10;
                PiHats_2(i,j) = meshgrids_ranges.p(i,j) * ( ...
                    (meshgrids_ranges.M(i,j) + 2/(gamma-1)) / (1 + 2/(gamma-1))...
                )^(2*gamma/(gamma-1));
            end            
            if length(PiCand) >=3
                warning("Three possible solutions detected using roots.")
            end
        elseif strcmpi('fsolve', backend_solver)
            if AR < 1
                [PiHats_2(i,j), ~, flagsPiHats_2(i,j)] = fsolve( ...
                @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                    Pi, ...
                    AR), 2, optimoptions('fsolve','Display','none'));
                if ~isreal(PiHats_2(i,j))
                    [PiHats_2(i,j), ~, flagsPiHats_2(i,j)] = fsolve( ...
                        @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                            Pi, ...
                            AR), 1, ...
                            optimoptions('fsolve','Display','none'));
                end
                if flagsPiHats_2(i,j) ~= 1
                    PiHats_2(i,j) = NaN;
                end
                [PiHats_sup(i,j), ~, flagsPiHats_sup(i,j)] = fsolve( ...
                    @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                        Pi, ...
                        AR), 0.4, ...
                        optimoptions('fsolve','Display','none'));
                if flagsPiHats_sup(i,j) ~= 1
                    PiHats_sup(i,j) = NaN;
                end
            elseif true % AR > 1
                [PiHats_2(i,j), ~, flagsPiHats_2(i,j)] = fsolve( ...
                @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                    Pi, ...
                    AR), 0.5, optimoptions('fsolve','Display','none'));
                [PiHats_sup(i,j), ~, flagsPiHats_sup(i,j)] = fsolve( ...
                @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                    Pi, ...
                    AR), 1.4, optimoptions('fsolve','Display','none'));
                if ~isreal(PiHats_2(i,j))
                    [PiHats_2(i,j), ~, flagsPiHats_2(i,j)] = fsolve( ...
                        @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                            Pi, ...
                            AR), real(PiHats_2(i,j)), ...
                            optimoptions('fsolve','Display','none'));
                end

                if flagsPiHats_2(i,j) ~= 1
                    PiHats_2(i,j) = NaN;
                end
                if flagsPiHats_sup(i,j) ~= 1
                    PiHats_sup(i,j) = NaN;
                end

            else
                [PiHats_2(i,j), ~, flagsPiHats_2(i,j)] = fzero( ...
                @(Pi) eqn_grid2target(meshgrids_ranges.data(i,j), ...
                    Pi, ...
                    AR), [1.5+1e-5,4]);
            end
        end
    end
end

% Map P-hat to M-hat
PiHat2MHat = @(Pi) sqrt( ...
        2/(gamma-1)*(1-Pi.^((gamma-1)/gamma)) ...
        ./(Pi.^((gamma-1)/gamma)-Pi.^((gamma+1)/gamma) ./ AR.^2));
MHats_2 = PiHat2MHat(PiHats_2);
MHats_2(flagsPiHats_2 == 10) = 1;
MHats_sups = PiHat2MHat(PiHats_sup);

% Create export
exports = struct( ...
    'Pi', PiHats_2, ...
    'M', MHats_2, ...
    'Pi_sups', PiHats_sup, ...
    'M_sups', MHats_sups,...
    'flags', flagsPiHats_2, ...
    'flags_sups', flagsPiHats_sup, ...
    'M_grid', meshgrids_ranges.M, ...
    'p_grid', meshgrids_ranges.p ...
);

%% Streamslice plot

function subs_plot()
% streamStartX = [linspace(0,2,10), 2+0*linspace(0,2,10)];
% streamStartY = [0*linspace(0,2,10), linspace(0,2,10)];
% sl_handle = streamline(meshgrids_ranges.M, meshgrids_ranges.p, ...
%     (MHats_2-meshgrids_ranges.M), ...
%     (PiHats_2-meshgrids_ranges.p), ...
%     streamStartX, streamStartY);
sl_handle = streamslice(meshgrids_ranges.M, meshgrids_ranges.p, ...
    (MHats_2-meshgrids_ranges.M), ...
    (PiHats_2-meshgrids_ranges.p), 2);
set(sl_handle,'color',[1 0.0 0.1]);
ylim([0 4])
xlim([0 2])

% Plot equilibrium contour
hold on
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    MHats_2-meshgrids_ranges.M, [0 0]);
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    PiHats_2-meshgrids_ranges.p, [0 0], 'LineStyle', '--', 'LineColor', [0 0 0]);
hold off

xlabel 'M_R'
ylabel '\Pi_R'
axis equal
end

function subs_clipped_plot()
% streamStartX = [linspace(0,2,10), 2+0*linspace(0,2,10)];
% streamStartY = [0*linspace(0,2,10), linspace(0,2,10)];
% sl_handle = streamline(meshgrids_ranges.M, meshgrids_ranges.p, ...
%     (MHats_2-meshgrids_ranges.M), ...
%     (PiHats_2-meshgrids_ranges.p), ...
%     streamStartX, streamStartY);

sonicPressure = @(p, M) ...
    p .* ...
    ((1 + (gamma-1)/2) ./ (1 + (gamma-1)/2 * M.^2)).^(-gamma / (gamma-1));
sonicPressures = sonicPressure(PiHats_2, MHats_2);

if AR < 1
    areaToSonic = 1/AR;
else
    areaToSonic = AR;
end

targetPressure = PiHats_2;
targetMach = MHats_2;

%% Compute choked-outlet target state
if AR <= 1
MSpecial = fzero( @(M) ...
    ((gamma+1)/2)^(-(gamma+1)/2/(gamma-1)) * ...
    (1 + (gamma-1)/2 * M^2 )^ ...
    ((gamma+1)/2/(gamma-1)) ./ M - ...
    areaToSonic, ...
    [1e-14,1-1e-14]);
    pSpecial = pressureMachFunction(gamma, MSpecial) / ...
           pressureMachFunction(gamma, 1);

    targetPressure (sonicPressures > 1) = pSpecial;
    targetMach (sonicPressures > 1) = MSpecial;
else
%     MExit = fzero( @(M) ...
%         ((gamma+1)/2)^(-(gamma+1)/2/(gamma-1)) * ...
%         (1 + (gamma-1)/2 * M^2 )^ ...
%         ((gamma+1)/2/(gamma-1)) ./ M - ...
%         areaToSonic, ...
%         [1e-14,1-1e-14]);
%     MSpecial = 1;
%     pSpecial = pressureMachFunction(gamma, 1) / ...
%            pressureMachFunction(gamma, MExit);

    for i = 1:size(targetPressure,1)
        for j = 1:size(targetPressure,2)
            if sonicPressures(i,j) > 1
                % Compute choked-inlet pressure from invariant
                targetMach(i,j) = 1;
                targetPressure(i,j) = meshgrids_ranges.p(i,j) * ( ...
                    (meshgrids_ranges.M(i,j) + 2/(gamma-1)) / (1 + 2/(gamma-1))...
                )^(2*gamma/(gamma-1));
            
                % Get root solver library
                [~, funcs] = bc_stability_utility();
                [pChokedMax, pChokedMin] = funcs.get_choked_pressure_range(AR, gamma);
                
                if targetPressure(i,j) > pChokedMax
                    targetPressure(i,j) = pChokedMax;
                end
            end
        end
    end
end

%% Seek lower pressure no-solution region
pFinder = [];
for j = 1:size(PiHats_2,2)
    pFinder(j) = meshgrids_ranges.p( ...
        find(arrayfun(@(x) ~isnan(x), PiHats_2(:,j)), 1, 'first'),...
        j);
end
% Replace lower pressure no-solution region with zero-flow target
for i = 1:size(meshgrids_ranges.p, 1)
    for j = 1:size(meshgrids_ranges.p, 2)
        if meshgrids_ranges.p(i,j) < pFinder(j)
            targetMach(i,j) = 0;
            targetPressure(i,j) = 1;
        end
    end
end

%% Create streamplot
sl_handle = streamslice(meshgrids_ranges.M, meshgrids_ranges.p, ...
    (targetMach-meshgrids_ranges.M), ...
    (targetPressure-meshgrids_ranges.p), 2);
set(sl_handle,'color',[1 0.0 0.1]);
ylim([0 4])
xlim([0 2])

% Plot equilibrium contour for both Pi solutions
hold on
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    MHats_2-meshgrids_ranges.M, [0 0]);
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    PiHats_2-meshgrids_ranges.p, [0 0], 'LineStyle', '--', 'LineColor', [0 0 0]);
hold off

if AR > 1
    % Plot locus of inlet-choked correction
    hold on
    contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
        targetPressure-meshgrids_ranges.p, [0 0], 'LineStyle', '--', 'LineColor', [0.5 0.5 0]);
    hold off
end

xlabel 'M_R'
ylabel '\Pi_R'
axis equal
end

function sups_plot()
% streamStartX = [linspace(0,2,10), 2+0*linspace(0,2,10)];
% streamStartY = [0*linspace(0,2,10), linspace(0,2,10)];
% sl_handle = streamline(meshgrids_ranges.M, meshgrids_ranges.p, ...
%     (MHats_2-meshgrids_ranges.M), ...
%     (PiHats_2-meshgrids_ranges.p), ...
%     streamStartX, streamStartY);
if AR < 1
    sl_handle = streamslice(meshgrids_ranges.M, meshgrids_ranges.p, ...
    (MHats_sups-meshgrids_ranges.M), ...
    (PiHats_sup-meshgrids_ranges.p), 2);
    set(sl_handle,'color',[0.1 0.0 1.0]);
    ylim([0 4])
    xlim([0 2])
else
    dM = MHats_sups-meshgrids_ranges.M;
    dp = PiHats_sup-meshgrids_ranges.p;
    % Partition data into physical (sufficient sonic pressure) and unphys.
    sonicPressure = @(p, M) ...
        p .* ...
        ((1 + (gamma-1)/2) ./ (1 + (gamma-1)/2 * M.^2)).^(-gamma / (gamma-1));
    sonicPressures = sonicPressure(PiHats_sup, MHats_sups);
    dM_phys = dM;
    dp_phys = dp;
    dM_phys(sonicPressures < 1) = NaN;
    dp_phys(sonicPressures < 1) = NaN;
    
    dM_unphys = dM;
    dp_unphys = dp;
    dM_unphys(sonicPressures >= 1) = NaN;
    dp_unphys(sonicPressures >= 1) = NaN;
    
    sl_handle = streamslice(meshgrids_ranges.M, meshgrids_ranges.p, ...
    dM_phys, ...
    dp_phys, 2);
    set(sl_handle,'color',[1 0.0 0.1]);
    
    hold on
    sl_handle = streamslice(meshgrids_ranges.M, meshgrids_ranges.p, ...
    dM_unphys, ...
    dp_unphys, 2);
    set(sl_handle,'color',[0.1 0.0 1]);
    hold off
    
    ylim([0 4])
    xlim([0 2])
end

% Plot equilibrium contour
hold on
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    MHats_sups-meshgrids_ranges.M, [0 0]);
contour(meshgrids_ranges.M, meshgrids_ranges.p, ...
    PiHats_sup-meshgrids_ranges.p, [0 0], 'LineStyle', '--', 'LineColor', [0 0 0]);
hold off

xlabel 'M_R'
ylabel '\Pi_R'
axis equal
end

exports.subs_plot = @()subs_plot();
exports.subs_clipped_plot = @()subs_clipped_plot();
exports.sups_plot = @()sups_plot();
exports.sups_clipped_plot = @()sups_plot();

end