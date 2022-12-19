% Add dependencies from repository root
addpath .\FlowRelations
addpath .\SBPSAT
addpath ..\sbplib

figure(7); clf;
agtools.plotOperatingDistribution(solution, metadata);