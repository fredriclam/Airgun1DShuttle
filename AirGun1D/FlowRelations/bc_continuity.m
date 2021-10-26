
%%
gamma = 1.4;
sonicContinuity = @(Pi, M, AR) ...
    (1 + (gamma-1)/2 * M.^2) .* Pi .^ ((gamma-1)/gamma) - ...
    (1 + (gamma-1)/2 * M.^2 .* AR.^(-2) .* Pi .^ ((gamma+1)/gamma));
sonicPressure = @(p, M) ...
    p .* ...
    ((1 + (gamma-1)/2) ./ (1 + (gamma-1)/2 * M.^2)).^(-gamma / (gamma-1));

% Principal square root
M2solve = @(Pi, AR) ...
    ( 2/(gamma-1) * (1 - Pi.^((gamma-1)/gamma)) ./ ...
          (Pi.^((gamma-1)/gamma) - Pi.^((gamma+1)/gamma)./ AR.^2)  );
Msolve = @(Pi, AR) ...
    sqrt( 2/(gamma-1) * (1 - Pi.^((gamma-1)/gamma)) ./ ...
          (Pi.^((gamma-1)/gamma) - Pi.^((gamma+1)/gamma)./ AR.^2)  );

%%

ranges = struct();
ranges.Pi = linspace(0.0, 2, 100*8);
ranges.AR = linspace(0.0, 2, 50*8);
meshgrids = struct();
[meshgrids.AR, meshgrids.Pi] = meshgrid(ranges.AR, ranges.Pi);

Mhats = nan(length(ranges.Pi), length(ranges.AR));
flags = nan(length(ranges.Pi), length(ranges.AR));
useRootSolver = false;

for i = 1:length(ranges.Pi)
    for j = 1:length(ranges.AR)
        if useRootSolver
            [Mhats(i,j), ~, flags(i,j)] = fsolve( ...
                @(M)sonicContinuity(ranges.Pi(i), M, ranges.AR(j)), .1, ...
                optimoptions('fsolve','Display','off','OptimalityTolerance', 1e-6));
            if flags(i,j) ~= 1
                [Mhats(i,j), ~, flags(i,j)] = fsolve( ...
                    @(M)sonicContinuity(ranges.Pi(i), M, ranges.AR(j)), 0.1, ...
                    optimoptions('fsolve','Display','off','OptimalityTolerance', 1e-6));
            end
        else
            Mhats(i,j) = M2solve(ranges.Pi(i), ranges.AR(j));
        end
    end
end
if useRootSolver
    Mhats(flags ~= 1) = nan;
else
    Mhats(~isreal(Mhats)) = nan;
    for i = 1:size(Mhats,1)
        for j = 1:size(Mhats,2)
            if Mhats(i,j) < 0
                Mhats(i,j) = nan;
            else
                Mhats(i,j) = sqrt(Mhats(i,j));
            end
        end
    end
end
sonicPressureRatios = sonicPressure(meshgrids.Pi, Mhats);

%% Draw (legacy)

if false
    figure(51);
    subplot(2,2,1);
    surf(ranges.AR, ranges.Pi, Mhats, 'LineStyle', 'none')
    xlabel 'AR'
    ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
    zlabel 'M'
    zlim([0, 1])
    drawnow
    
    subplot(2,2,2);
    surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
    xlabel 'AR'
    ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
    zlabel ('${p}^* / p_b$','Interpreter', 'latex')
    zlim([0, 10])
    
    Mhats_subsbranch = Mhats;
    Mhats_subsbranch(sonicPressureRatios>1) = nan;
    Mhats_subsbranch(Mhats_subsbranch > 1 | Mhats_subsbranch < 1e-3) = nan;
    
    subplot(2,2,3);
    surf(meshgrids.AR, meshgrids.Pi, Mhats_subsbranch, 'LineStyle', 'none')
    % hold on
    % surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
    xlabel 'AR'
    ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
    zlabel 'M'
    zlim([0, 2])
    title 'M filtered to p^* < p_b, M < 1'
    view([0 0 1])
    colorbar
    caxis([0, 1])
    xlim([0.5, 3])
    ylim([0, 2])
    
    subplot(2,2,4);
    mhats_supsbranch = Mhats;
    mhats_supsbranch(sonicPressureRatios<1) = nan;
    
    surf(meshgrids.AR, meshgrids.Pi, mhats_supsbranch, 'LineStyle', 'none')
    % hold on
    % surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
    xlabel 'AR'
    ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
    zlabel 'M'
    zlim([0, 2])
    title 'M > 1?'
    view([0 0 1])
    colorbar
    caxis([0, 2])
    xlim([0.5, 3])
    ylim([0, 2])
end

%% Export graph
figure(52); clf;

clim = [0, 1];
Mhats_clipped = Mhats;
Mhats_clipped(Mhats_clipped > clim(2)) = clim(2);
contourf(meshgrids.AR, meshgrids.Pi, Mhats_clipped, ...
    'LineStyle', 'none', 'LevelStep', 0.05);
% hold on
% surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
xlabel 'AR'
ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
title 'M'
cb = colorbar();
cb.Label.String = 'M';

caxis(clim);


xlim([0, 2])
ylim([0, 2])

%%
hold on
contour(meshgrids.AR, meshgrids.Pi, Mhats, [1, 1], ...
        'linecolor', [0 1 0], 'linewidth', 1.5)
contour(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, [1, 1], ...
        'linecolor', [0.8 0.5 0.5], 'linewidth', 1.5)
plot(linspace(0.,2.,100), linspace(0.,2.,100).^1.4, ...
    'Color', [0.8, 0.1, 0.2], 'linewidth', 1.5)
hold off

legend( ...
    ["$\hat{M}$ surface", ...
    "$\hat{M} = 1$", ...
    "$\hat{p}^* / p_\mathrm{b}$", ...
    "$\mathbf{AR}^\gamma$"], ...
    'Interpreter', 'latex')
axis equal

%%
figure(53) ;clf;

clim = [0, 1];
sonicPressureRatios_clipped = sonicPressureRatios;
sonicPressureRatios_clipped(sonicPressureRatios > clim(2)) = clim(2);
contourf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios_clipped, ...
    'LineStyle', 'none', 'LevelStep', 0.05);
% hold on
% surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
xlabel 'AR'
ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
title ('$\hat{p^*} / p_b$','Interpreter', 'latex')
cb = colorbar();
cb.Label.String = '$\hat{p^*} / p_b$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

caxis(clim);

xlim([0, 2])
ylim([0, 2])

hold on
contour(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, [1, 1], ...
        'linecolor', [0.8 0.5 0.5], 'linewidth', 1.5)
hold off

%% Physical selection
figure(54); clf;
Mphysical = Mhats;
Mphysical(sonicPressureRatios > 1) = nan;
Mphysical(Mhats > 1) = nan;

clim = [0, 1];
Mphysical_clipped = Mphysical;
Mphysical_clipped(Mphysical > clim(2)) = clim(2);
contourf(meshgrids.AR, meshgrids.Pi, Mphysical_clipped, ...
    'LineStyle', 'none', 'LevelStep', 0.05);
% hold on
% surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
xlabel 'AR'
ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
title ('$M$, physical regions','Interpreter', 'latex')
cb = colorbar();
cb.Label.String = '$\hat{M}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

caxis(clim);

xlim([0, 2])
ylim([0, 2])

sonicPressureRatios_copy = sonicPressureRatios;
sonicPressureRatios_copy(meshgrids.Pi < 1) = nan;
Mhats_copy = Mhats;
Mhats_copy(meshgrids.Pi < 1) = nan;

hold on
contour(meshgrids.AR, meshgrids.Pi, sonicPressureRatios_copy, [1, 1], ...
        'linecolor', [0, 0, 0], 'linewidth', 1.5)
contour(meshgrids.AR, meshgrids.Pi, Mhats_copy, [1, 1], ...
        'linecolor', [0, 0, 0], 'linewidth', 1.5)
hold off


%% Use model results to compare
figure(55); clf;
Mphysical = Mhats;
Mphysical(sonicPressureRatios > 1) = nan;
Mphysical(Mhats > 1) = nan;

clim = [0, 1];
Mphysical_clipped = Mphysical;
Mphysical_clipped(Mphysical > clim(2)) = clim(2);
contourf(meshgrids.AR, meshgrids.Pi, Mphysical_clipped, ...
    'LineStyle', 'none', 'LevelStep', 0.05);
% hold on
% surf(meshgrids.AR, meshgrids.Pi, sonicPressureRatios, 'LineStyle', 'none')
xlabel 'AR'
ylabel ('$\hat{p} / p_b$','Interpreter', 'latex')
title ('$M$, physical regions','Interpreter', 'latex')
cb = colorbar();
cb.Label.String = '$\hat{M}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

caxis(clim);

xlim([0, 2])
ylim([0, 2])

sonicPressureRatios_copy = sonicPressureRatios;
sonicPressureRatios_copy(meshgrids.Pi < 1) = nan;
Mhats_copy = Mhats;
Mhats_copy(meshgrids.Pi < 1) = nan;

hold on
contour(meshgrids.AR, meshgrids.Pi, sonicPressureRatios_copy, [1, 1], ...
        'linecolor', [0, 0, 0], 'linewidth', 1.5)
contour(meshgrids.AR, meshgrids.Pi, Mhats_copy, [1, 1], ...
        'linecolor', [0, 0, 0], 'linewidth', 1.5)
hold off

pS =[fullState.portStates];
bS = [fullState.bubbleStates];
hold on
plot( [pS.APortExposed] ./ ...
          (metadata.paramAirgun.airgunCrossSecAreaSqInch * 0.0254^2), ...
      [pS.pPort] ./ [bS.p], ...
      'm-', 'LineWidth', 1.5);
hold off