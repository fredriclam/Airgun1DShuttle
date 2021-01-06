% Plotting scripts for extra experimental comparisons
% Test case 83
% Run after test_launch_script_v1

%% Right-side pressure comparison
figure(1021); clf;

subplot(1,3,1);
plot(sol.x, p(end,:));
xlim([0, 0.4]);
ylim([0, 8e6]);
title 'Pressure at right-end PDE domain (fixed time)'

subplot(1,3,2);
plot(sol.x, p2(end,:));
xlim([0, 0.4]);
ylim([0, 8e6]);
title 'Pressure at right-end PDE domain (w/ shuttle)'

subplot(1,3,3);
% Unscaled
% plot(HiTestData(24).iNetTimeAxisP, HiTestData(24).iNetCh13Data)
% Scaled
startTimeIndex = 18400;
endTimeIndex = startTimeIndex+2000;
plot(HiTestData(24).iNetTimeAxisP(startTimeIndex:endTimeIndex) -  ...
     HiTestData(24).iNetTimeAxisP(startTimeIndex), ...
     6894.76 * (1050 + HiTestData(24).iNetCh13Data(startTimeIndex:endTimeIndex)));
ylim([0, 8e6]);
xlabel([0, 0.4]);
title 'Experimental pressure (id:24, Ch:13) just inside port'
ylabel 'Pa (shifted vert and horizontally) for gauge pressure comp'
xlabel 't [s]'

%% Closed-end
%% Right-side pressure comparison
figure(1022); clf;

subplot(1,3,1);
plot(sol.x, p(1,:));
xlim([0, 0.4]);
ylim([0, 8e6]);
title 'Pressure at closed-end PDE domain (fixed time)'

subplot(1,3,2);
plot(sol.x, p2(1,:));
xlim([0, 0.4]);
ylim([0, 8e6]);
title 'Pressure at closed-end PDE domain (w/ shuttle)'

subplot(1,3,3);
% Unscaled
% plot(HiTestData(24).iNetTimeAxisP, HiTestData(24).iNetCh13Data)
% Scaled
startTimeIndex = 18400;
endTimeIndex = startTimeIndex+2000;
plot(HiTestData(24).iNetTimeAxisP(startTimeIndex:endTimeIndex) -  ...
     HiTestData(24).iNetTimeAxisP(startTimeIndex), ...
     6894.76 * (1000 + HiTestData(24).iNetCh16Data(startTimeIndex:endTimeIndex)));
ylim([0, 8e6]);
xlabel([0, 0.4]);
title 'Experimental pressure (id:24, Ch:16) just inside port'
ylabel 'Pa (shifted vert +1000psi and horizontally) for gauge pressure comp'
xlabel 't [s]'

%% SNAP data
figure(1023); clf;

% subplot(1,3,1);
% plot(sol.x, p(1,:));
% xlim([0, 0.4]);
% ylim([0, 8e6]);
% title 'Pressure at closed-end PDE domain (fixed time)'
% 
% subplot(1,3,2);
% plot(sol.x, p2(1,:));
% xlim([0, 0.4]);
% ylim([0, 8e6]);
% title 'Pressure at closed-end PDE domain (w/ shuttle)'



% Do bubble
figPos = get(gcf,'Position');
% set(gcf,'Position',[figPos(1) 0 600 900]);

subplot(1,2,1);

% bubble volume
V = 4/3*pi*R.^3;
V2 = 4/3*pi*R2.^3;


% acoustic pressure
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r, 'LineWidth', 1);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
hold on;
plot((tInterp2-r/c_inf)*1000, pPres2*1e-5*r, 'k', 'LineWidth', 1);
legend({'Instant open 10 ms','Shuttle'},'location','best')
xlim([0, 400])


subplot(1,2,2);
% Unscaled
% plot(HiTestData(24).iNetTimeAxisP, HiTestData(24).iNetCh13Data)
% Scaled
startTimeIndex = HiTestData(24).t0indexGreen-1000;
endTimeIndex = startTimeIndex+32000*0.4;
dt = 1/HiTestData(24).fsSNAP;
plot(dt*(0:endTimeIndex-startTimeIndex), ...
     HiTestData(24).dataGreen(startTimeIndex:endTimeIndex));
xlim([0, 0.4])
xlabel([0, 0.4]);
title 'Experimental pressure signal (id:24), SNAP'
ylabel 'SNAP data, voltage [V]'
xlabel 't [s]'