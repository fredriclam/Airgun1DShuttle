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
M_RHistory = [eDS.M_R];

pS = [gridFullStates.portStates];
portArea = [pS.APortExposed];
csArea = metadata.discretization.physConst.crossSectionalArea;

bS = [gridFullStates.bubbleStates];
pBubbleHistory = [bS.p];

pRatio = pSonic_RHistory ./ pBubbleHistory ;
ARatio = portArea / csArea;

caseKeyHistory = cellfun(@(k) caseKey2Num(k), {pS.caseKey});

figure(1); clf;
subplot(1,2,1);

if false
    % Plain color plot
    plot(pRatio(1), ARatio(1), 'k.', 'MarkerSize', 24)
    text(pRatio(1),0.05,'Start')
    hold on
    plot(pRatio, ARatio, 'k', 'LineWidth', 1)
    hold off
else
    % Color coded plot
    
    
    
    % Find all case key switches
    caseKeySwitchIndices = [0, find(...
        caseKeyHistory(2:end) ~= caseKeyHistory(1:end-1),...
        length(caseKeyHistory)-1)];
    
    colorMap = {'k', 'g', 'b', 'm', 'r', 'c'};
    % Dummy lines for legend
    for i = 1:5
        plot(pRatio(1), ARatio(1), [colorMap{i}, '-'])
        hold on
    end
    
    for i = 1:length(caseKeySwitchIndices)-1
        % Include the right boundary element too for continuity
        plotRange = caseKeySwitchIndices(i)+1:caseKeySwitchIndices(i+1)+1;
        plot(pRatio(plotRange), ARatio(plotRange), ...
            [colorMap{1+caseKeyHistory(plotRange(1))}], 'LineWidth', 1);
        hold on
    end
    
    plot(pRatio(1), ARatio(1), 'k.', 'MarkerSize', 24)
    text(pRatio(1),0.05,'Start')
    hold off
    
    legendLabels = {'Closed', ...
        'Subsonic', ...
        'Port choked', ...
        'Chamber choked', ...
        'Chamber choked*'};
    legend(legendLabels, 'Interpreter', 'latex')
end

xlabel ('$p_\mathrm{R}^*/p_\mathrm{b}$', 'Interpreter', 'latex', 'FontSize', 18)
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

subplot(1,2,2);
plot(M_RHistory(1), ARatio(1), 'k.', 'MarkerSize', 24)
text(M_RHistory(1),0.05,'Start')
hold on
plot(M_RHistory, ARatio, 'k', 'LineWidth', 1)
hold off
xlabel ('$M_\mathrm{R}$', 'Interpreter', 'latex', 'FontSize', 18)
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

set(gcf, 'position', ...
    [-1714         369        1420         573]);

%% Figure 2-4: Wing plots
figure(2); clf;
plot(pRatio, solution.soln.x, 'k', 'LineWidth', 1)
xlabel ('$p_\mathrm{R}^*/p_\mathrm{b}$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel ('$t [s]$', 'Interpreter', 'latex', 'FontSize', 18)
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
set(gcf, 'position', ...
    [-1613          56         616         224]);

figure(3); clf;
plot(solution.soln.x, ARatio, 'k', 'LineWidth', 1)
xlabel ('$t [s]$', 'Interpreter', 'latex', 'FontSize', 18)
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
set(gcf, 'position', ...
    [-1903         376         306         565]);

figure(4); clf;
plot(M_RHistory, solution.soln.x, 'k', 'LineWidth', 1)
xlabel ('$M_\mathrm{R}$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel ('$t [s]$', 'Interpreter', 'latex', 'FontSize', 18)
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
set(gcf, 'position', ...
    [ -986    54   623   224]);

%% Figure 5: case keys
figure(5); clf;
caseKeyHistory = cellfun(@(k) caseKey2Num(k), {pS.caseKey});
plot(solution.soln.x, caseKeyHistory, 'k', 'LineWidth', 1)

xlabel ('$t [s]$', 'Interpreter', 'latex', 'FontSize', 18)
ylim([-0.5, 4.5])
set(gca, ...
    'FontSize', 14, ...
    'TickLabelInterpreter', 'latex', ...
    'YMinorTick', 'on', ...
    'XGrid', 'on', ...
    'XMinorGrid', 'on' ...
);

labels = {
    'Port closed', ...
    'Subsonic', ...
    'Port choked', ...
    'Chamber choked', ...
    'Chamber choked*'
};
set(gca, 'YTick', 0:4, 'YTickLabel', labels)

%% Explicit [deprecated]
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

function caseNum = caseKey2Num(caseKey)
    if strcmpi('portClosed', caseKey)
        caseNum = 0;
    elseif strcmpi('subsonic', caseKey)
        caseNum = 1;
    elseif strcmpi('portChoked', caseKey)
        caseNum = 2;
    elseif strcmpi('chamberChokedNatural', caseKey)
        caseNum = 3;
    elseif strcmpi('chamberChokedForced', caseKey)
        caseNum = 4;
    else
        caseNum = NaN;
    end
end