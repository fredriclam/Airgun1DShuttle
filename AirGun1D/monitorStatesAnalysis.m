% Clean up monitor data
if ~isfield(monitorStates, 'cleaned')
    t_raw = [monitorStates.data.t];
    ind_garbage = find(t_raw == max(t_raw), 1, 'first');
    monitorStates.data = monitorStates.data(1:ind_garbage);
    monitorStates.cleaned = true;
end

% Local processing
caseNums = [monitorStates.data.caseNum];

% Test plots
% plot([monitorStates.data.t], [monitorStates.data.bubbleMass]);
% plot([monitorStates.data.t], [monitorStates.data.pMid]);

% Filtering Euler domain states
rho2 = q2(1:3:end,:);
u2 =  q2(2:3:end,:) ./ rho2;
c_v = 718;
T2 =  (q2(3:3:end,:) ./ rho2 - 0.5 * u2.^2) / c_v;

% More filtering
% Hack: since the monitor data comes from each call from ODE45, some calls
% are discarded due to tolerance. Create index into monitorStates.data
% that preserves only the increasing t
ascendingFilter = [1];
tMaxCurrent = monitorStates.data(1).t;
for i = 2:length([monitorStates.data.t])
    if monitorStates.data(i).t <= tMaxCurrent
    
    else
        tMaxCurrent = monitorStates.data(i).t;
        ascendingFilter = [ascendingFilter, i];
    end
end

% Boundary condition type
% 0: closed
% 1: subsonic
% 2: sonic at port
% 3: sonic in firing chamber
monitorStatet = [monitorStates.data.t];
caseNumsODE = interp1(monitorStatet(ascendingFilter), ...
                      caseNums(ascendingFilter), t, 'nearest');

%% Plot 1: Exit condition history
plot(1e3*t, u2(end, :), 'k', 'LineWidth', 1);
hold on
plot(1e3*t(caseNumsODE == 0), u2(end, caseNumsODE == 0), '.m', 'MarkerSize', 10);
plot(1e3*t(caseNumsODE == 1), u2(end, caseNumsODE == 1), '.b', 'MarkerSize', 10);
plot(1e3*t(caseNumsODE == 2), u2(end, caseNumsODE == 2), '.g', 'MarkerSize', 10);
plot(1e3*t(caseNumsODE == 3), u2(end, caseNumsODE == 3), '.r', 'MarkerSize', 10);


hold off

xlabel 't [ms]'
ylabel 'Chamber exit velocity [m/s]'
legend({'all', 'closed', 'subsonic', 'sonic: port', 'sonic: chamber'})