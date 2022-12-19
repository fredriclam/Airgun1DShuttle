% Plot energy distribution in model
% (Hint: delete legend to show plot in regular aspect ratio)
figure(13); clf;
[~, exports] = agtools.plotEnergyDistro( ...
    solution_reference, metadata_reference);
set(gcf,'position',[360   103   716   841]);