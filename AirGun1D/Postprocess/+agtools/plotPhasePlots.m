% Plot 
function [ARatio, pRatio, t_] = plotPhasePlots(solution, metadata)

q_R = solution.q(end-2:end,:);
p_R = metadata.discretization.schm.p(solution.q(end-2:end,:));
c_R = metadata.discretization.schm.c(q_R);
u_R = q_R(2,:) ./ q_R(1,:);
M_R = u_R ./ c_R;

gamma = metadata.discretization.physConst.gamma;
% Fraction of stagnation pressure
fracStagation = @(M) ( 1 + (gamma-1)/2 * M.^2) .^ (-gamma / (gamma-1));
pSonic_R = fracStagation(1) ./ fracStagation(M_R) .* p_R;

t_ = solution.soln.x;

V_ = 4/3*pi*(solution.bubble(1,:)).^3;
E_ = solution.bubble(4,:);
pBubble_ = (gamma - 1) .* E_ ./ V_;

% plot(t_, pSonic_R ./ pBubble_);

portLength_ = metadata.discretization.physConst.airgunPortLength;
portLead_ = metadata.discretization.physConst.portLead;
APortTotal_ = metadata.discretization.physConst.APortTotal;
crossSectionalArea_ = ...
    metadata.discretization.physConst.crossSectionalArea;
APortExposed = (solution.shuttle(1,:) - portLead_) ...
    / (0.0254*portLength_) * APortTotal_;
% Clamp to [0, APortTotal]
APortExposed( APortExposed < 0 ) = 0;
APortExposed( APortExposed > APortTotal_ ) = APortTotal_;

% plot(t_, APortExposed ./ crossSectionalArea_);

% Plot trace
ARatio = APortExposed ./ crossSectionalArea_;
pRatio = pSonic_R ./ pBubble_;

tL = tiledlayout(2,2);
nexttile;
plot(APortExposed ./ crossSectionalArea_, pSonic_R ./ pBubble_);

xlabel ('$A_\mathrm{port} / A_\mathrm{cs}$', ...
        'Interpreter', 'latex', 'FontSize', 14)    
ylabel ('$p^*(x = 0) / p_\mathrm{b}$', ...
        'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'YScale', 'log')
grid on

nexttile;
plot(solution.bubbleContinuationState(1,:), ...
     solution.bubbleContinuationState(2,:))
hold on
plot(solution.bubble(1,end), ...
     solution.bubble(2,end), 'k.', 'MarkerSize', 12);
index = find(solution.bubbleContinuationTime > 1, 1, 'first');
tIndex = solution.bubbleContinuationTime(index);
plot(solution.bubbleContinuationState(1,index), ...
     solution.bubbleContinuationState(2,index), ...
     'k^', 'MarkerSize', 4);
legend(["Phase portrait", ...
        "t = 300 ms (port closed)",...
        "t = " + 1e3*tIndex + " ms"])
    
xlabel ('${R}$ (m)', 'Interpreter', 'latex', 'FontSize', 14)    
ylabel ('$\dot{R}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14)

nexttile;
VFn = @(R) 4/3*pi*R.^3;
VDotFn = @(R, Rdot) 4*pi*R.^2 .* Rdot;

% Compute V, Vdot for continued bubble model
VContinued = VFn(solution.bubbleContinuationState(1,:));
VDotContinued = VDotFn(solution.bubbleContinuationState(1,:), ...
                       solution.bubbleContinuationState(2,:));

plot(VContinued, VDotContinued);
xlabel ('${V}$ (m${}^3$)', 'Interpreter', 'latex', 'FontSize', 14)    
ylabel ('$\dot{V}$ (m${}^3$/s)', 'Interpreter', 'latex', 'FontSize', 14)
hold on
plot(VFn(solution.bubble(1,end)), ...
     VDotFn(solution.bubble(1,end), solution.bubble(2,end)), ...
     'k.', 'MarkerSize', 12);
plot(VContinued(index), ...
     VDotContinued(index), ...
     'k^', 'MarkerSize', 4);
legend(["Phase portrait", ...
        "t = 300 ms (port closed)",...
        "t = " + 1e3*tIndex + " ms"], 'location', 'eastoutside')

nexttile
% plot3(VContinued, VDotContinued, solution.bubbleContinuationState(4,:), ...
%     'LineWidth', 1.5);
patch([VContinued, NaN], ...
      [VDotContinued, NaN], ...
      [solution.bubbleContinuationState(4,:), NaN]/1e6, ...
      [solution.bubbleContinuationState(4,:), NaN]/1e6, ...
      'EdgeColor','interp','FaceColor','none','LineWidth',1.5)
view([1 1 1])
cb = colorbar;
cb.Label.String = 'E (MJ)';
xlabel ('${V}$ (m${}^3$)', 'Interpreter', 'latex', 'FontSize', 14)    
ylabel ('$\dot{V}$ (m${}^3$/s)', 'Interpreter', 'latex', 'FontSize', 14)
zlabel ('$E$ (MJ)', 'Interpreter', 'latex', 'FontSize', 14)
grid on
grid minor
colormap hsv
end