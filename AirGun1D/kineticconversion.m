% Kinetic energy conversion modeling
linearsys.totalp = 1e-6*postprocessStates(4).bubblePressure;
linearsys.t = 1000*postprocessStates(4).t;

k = linearsys.totalp(1);
p = 0;

for i = 2:length(linearsys.t)
    dt = linearsys.t(i) - linearsys.t(i-1);
    transferRate = 0.4;
    p(i) = p(i-1) + transferRate * (linearsys.totalp(i-1) - p(i-1)) * dt;
end

figure(10); clf;
plot(linearsys.t, linearsys.totalp, 'LineWidth', 1.5)
xlim([0, 20])
hold on
plot(linearsys.t, p, 'LineWidth', 1.5)