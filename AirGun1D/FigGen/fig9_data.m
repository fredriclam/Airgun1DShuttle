figure(9); clf;
subplot(3,1,1);
[t, p, pData] = agtools.plotFiringChamber_Wall(fullState, 'p', HiTestData);
xlabel('{\it{t}} (ms)', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('{\it{p}} (MPa)', 'FontSize', 14, 'Interpreter', 'tex');
subplot(3,1,2);
[t, T, TData] = agtools.plotFiringChamber_Wall(fullState, 'T', HiTestData);
xlabel('{\it{t}} (ms)', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('{\it{T}} (K)', 'FontSize', 14, 'Interpreter', 'tex');

subplot(3,1,3);
get_exp_entropy = @(p, T) ...
    T.^metadata_reference.discretization.physConst.gamma ./ ...
    p.^(metadata_reference.discretization.physConst.gamma-1);
es0 = get_exp_entropy(p(1,1), T(1,1));
get_entropy = @(p, T) log(get_exp_entropy(p,T) ./ es0);

plot(1e3*t, get_entropy(p, T) / es0, 'k--', 'LineWidth', 1);
hold on
plot(1e3*t, get_entropy(1e6*pchip(pData.t, pData.y, 1e3*t), ...
                        pchip(TData.t, TData.y, 1e3*t)), ...
     'b-', 'LineWidth', 1);
hold off

% Manual format
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('{\it{t}} (ms)', 'FontSize', 14);
ylabel("$\Delta s / c_v$", ...
       'Interpreter', 'latex', 'FontSize', 11);
legend({'Model', 'Data'}, ...
      'Interpreter', 'latex', 'FontSize', 13);
  
%%
for i = 1:3
    subplot(3,1,i);
    set(legend, 'Interpreter', 'tex');
    set(gca, 'TickLabelInterpreter', 'tex', 'FontSize', 14)
end