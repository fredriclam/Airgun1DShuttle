q1History = solution.q(1:3:end,:);
q2History = solution.q(2:3:end,:);
q3History = solution.q(3:3:end,:);


tHistory = solution.soln.x;

% (tHistory, q1History)
% 
% metadata.discretization.schm.p()

% Initialize scalar values at eval points
rhoHistory = nan(size(q1History));
pHistory = nan(size(q1History));
uHistory = nan(size(q1History));
cHistory = nan(size(q1History));
MHistory = nan(size(q1History));
THistory = nan(size(q1History));

clear gridFullStates;
%% Full state function
% for i = 1:size(q1History,1)
for j = 1:size(q1History,2)
    bubble = solution.bubble(:,j);
    shuttle = solution.shuttle(:,j);
    t = solution.soln.x(j);
    gridFullStates(j) = metadata.discretization.fullState(...
        solution.q(:,j), ...
        t, ...
        bubble, ...
        shuttle, ...
        false ...
    );
end 
% end

%% Figure 1: Phase plot
eDS = [gridFullStates.eulerDomainStates];
pSonic_RHistory = [eDS.pSonic_R];
p_RHistory = [eDS.p_R];

pS = [gridFullStates.portStates];
portArea = [pS.APortExposed];
csArea = metadata.discretization.physConst.crossSectionalArea;

%
pRatio = p_RHistory ./ pSonic_RHistory;
ARatio = portArea / csArea
plot(pRatio(1), ARatio(1), 'k.', 'MarkerSize', 24)
hold on
plot(pRatio, ARatio, 'k', 'LineWidth', 1)
hold off
xlabel ('$p_\mathrm{R}/p_\mathrm{R}^*$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel ('$A_\mathrm{port}/A_\mathrm{cs}$', 'Interpreter', 'latex', 'FontSize', 18)
set(gca, ...
    'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
	'XGrid', 'on', ...
    'XMinorGrid', 'on', ...
    'YGrid', 'on', ...
    'YMinorGrid', 'on' ...
);



%% Explicit
% Extract states manually
% for i = 1:size(q1History,1)
%     % Extract q for this grid point x, for all t
%     qHistoryFixedSpace = [q1History(i,:); 
%                           q2History(i,:);
%                           q3History(i,:)];
%     % Compute useful states
%     rhoHistory(i,:) = q1History(i,:);
%     pHistory(i,:) = metadata.discretization.schm.p(qHistoryFixedSpace);
%     uHistory(i,:) = q2History(i,:) ./ q1History(i,:);
%     cHistory(i,:) = metadata.discretization.schm.c(qHistoryFixedSpace);
%     MHistory(i,:) = uHistory(i,:) ./ cHistory(i,:);
%     THistory(i,:) = pHistory(i,:) ./ ...
%         (metadata.discretization.physConst.Q * rhoHistory(i,:));
% end