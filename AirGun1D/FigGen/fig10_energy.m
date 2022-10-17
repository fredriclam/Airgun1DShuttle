% Energy distribution
% (Delete legend to show plot in regular aspect ratio)

figure(10); clf;
[~, exports] = agtools.plotEnergyDistro(solution_reference, metadata_reference);

set(gcf,'position',[360   103   716   841]);