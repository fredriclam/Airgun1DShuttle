% Shuttle coupled to acoustic equation with
% u(0, t) boundary condition at xi = 0, and a wall at xi = -L, i.e.,
% u(-L, t) = 0; dP/dx(-L, t) = 0.
% DDE:
%   m d^2 xi / dt^2 + (rho0 c) * (dxi/dt + dxi/dt(t - 2L/c)) + A_R * pR(xi)
%      = A_L * p_{L,0}

%% Params
clear;
m = 63 * .454;                % [kg]
% L = 0.6;                      % [m]
p0 = 2000 * 0.068046e5;       % [Pa]
T0 = 300;                     % [K]
gamma = 1.4;                  % [-]
R = 287;                      % [J / kg K]
c = sqrt(gamma * R * T0);     % [m/s]
rho0 = p0 / (R * T0);         % [kg/m^3]

cushLength = 3.0 * 0.0254;    % [m]
accelLength = 0.0 * 0.0254;   % [m]
pRFactor = @(xi) (xi <= accelLength) .* (1) + ...
    (xi > accelLength) .* ...
    (cushLength ./ (cushLength - (xi - accelLength))) .^ gamma; % [-]
pR = @(xi) p0 * pRFactor(xi); % [Pa]

A_L = 10 * 16 * (0.0254^2);    % [m^2] -- from test_launch_script test case 1
A_R = 0.0506/0.0586 * A_L;    % [m^2]

L = c * pi / sqrt((A_R * p0 * gamma / cushLength)/m); % Resonant length [m]

cushMass = p0 / (R * T0) * A_R * cushLength; % [kg]
c_v = 1/(gamma-1)*R;                         % [J / kg K]
energyDelta = @(xi) (pRFactor(xi).^((gamma-1)/gamma) - 1) * ...
    cushMass * c_v * T0;

%% Compute equilibrium position
forceFactor = ((A_L * p0) / (A_R * p0))^(1/gamma);
xi_eq = (forceFactor-1)/forceFactor * cushLength;

%% Solve and plot
figure(1); clf;
tFinal = 0.100;

sweepFactorVector = [0.1, linspace(2/3,3/2,1), 10];
collection = [];
for i = 1:length(sweepFactorVector)
    [tVector, xVector] = DDESolve(tFinal, true, true, ...
        m, sweepFactorVector(i) * L, p0, T0, gamma, R, c, rho0, ...
        cushLength, accelLength, pR, A_R, A_L);
    collection(1).tVector = tVector;
    collection(1).xVector = xVector;
    
    subplot(1,3,1);
    plot(tVector, xVector(1,:));
    if i == 1
        hold on
    end
    title('L ~ L_{res}')
    xlabel('time [s]')
    ylabel('shuttle position [m]')
    
    subplot(1,3,2);
    plot(tVector, energyDelta(xVector(1,:)) + ...
        0.5*m*xVector(2,:).^2);
    if i == 1
        hold on
    end
    title ('Shuttle + cushion energy')
    
    
    % hold on;
    % subplot(2,3,2);
    % plot(tVector1, xVector1(2,:));
    % hold on;
    % title('Velocity [m/s] vs. time [s]')
    subplot(1,3,3);
    dt = tVector(2)-tVector(1);
    [f, X] = outil.spect(xVector(1,:), dt);
    loglog(f, X);
    if i == 1
        hold on
    end
    drawnow;
end

% % L < Resonant length solve
% [tVector2, xVector2] = DDESolve(tFinal, true, true, ...
%     m, 0.1*L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, A_R, A_L);
% % Primary rarefaction plot
% subplot(2,3,2);
% plot(tVector2, xVector2(1,:));
% title('L ~ 0.1L_{res}')
% % subplot(1,2,2);
% % plot(tVector2, xVector2(2,:));
% subplot(2,3,5);
% dt = tVector2(2)-tVector2(1);
% [f2, X2] = outil.spect(xVector2(1,:), dt);
% loglog(f2, X2);
% 
% % L > Resonant length solve
% [tVector3, xVector3] = DDESolve(tFinal, true, true, ...
%     m, L*10, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, A_R, A_L);
% % Damped plot
% subplot(2,3,3);
% plot(tVector3, xVector3(1,:));
% title('L ~ 10L_{res}')
% % subplot(1,2,2);
% % plot(tVector3, xVector3(2,:));
% subplot(2,3,6);
% dt = tVector3(2)-tVector3(1);
% [f3, X3] = outil.spect(xVector3(1,:), dt);
% loglog(f3, X3);
% 
% % % Guidelines
% % subplot(1,2,1);
% % % Equilibrium by force calulation
% % plot([0, tFinal], accelLength + xi_eq * [1, 1], 'k--');
% % % Shuttle round-trip times
% % yyy = ylim;
% % plot(0.005*[1 1], [yyy(1)-1000, yyy(2)+1000], 'k--');
% % plot(0.008*[1 1], [yyy(1)-1000, yyy(2)+1000], 'k--');
% % ylim(yyy);
% % 
% % subplot(1,2,2);
% % % Shuttle round-trip times
% % yyy = ylim;
% % plot(0.005*[1 1], [yyy(1)-1000, yyy(2)+1000], 'k--');
% % plot(0.008*[1 1], [yyy(1)-1000, yyy(2)+1000], 'k--');
% % ylim(yyy);
% % 
% % legend({'L_res', 'L_res/2', '2L_res'}, ...
% %     'Location', 'northeast')

%% Figure 2: frequency superposition
% figure(2);
% clf
% loglog(f1,X1); hold on; loglog(f2,X2); loglog(f3,X3);

%% Function: DDE solve with fixed timestep
function [tVector, xVector] = DDESolve(tFinal, globalDamp, lagDamp, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, A_R, A_L)

tFinal = 0.100;
N = 1000;
tVector = linspace(0, tFinal, N+1);
xVector = nan(2, length(tVector));
dt = tFinal/N;

% Set initial conditions
xVector(:,1) = [0; 0];

for i = 1:N
    t = tVector(i);
    x = xVector(:, i);
    assert(~any(isnan(x)) && all(isreal(x)))
    
    
    
    % Sum contributions over domain of influence:
    damping = x(2) + ...
        imageContrib(t, L, c, lagDamp, ...
        tVector, xVector);
    
    % Build RHS at time t
    dx = [x(2); ...
          1/m * (...
              A_L * p0 - A_R * pR(x(1)) ...
              - globalDamp * rho0*c * A_L * damping)
         ] * dt;
    
    xGuess = x + dx / 2;
    
    % Build approximate RHS at time t + dt/2; use same xVector (stability
    % condition needs to be respected)
    damping = xGuess(2) + ...
        imageContrib(t + dt/2, L, c, lagDamp, ...
        tVector, xVector);
    
    dx = [xGuess(2); ...
          1/m * (...
              A_L * p0 - A_R * pR(xGuess(1)) ...
              - globalDamp * rho0*c * A_L * damping)
         ] * dt;
    % Get next x from Forward Euler operator
    xVector(:, i+1) = x + dx;
end
end

%% Function: samples a data history at query time
% Input:
%   t -- Domain vector
%   data -- Image vector
%   tq -- Query time
function y = sampleHistory(t, data, tq)
    if tq < 0
        y = 0;
    else
        y = interp1(t, data, tq);
    end
    assert(~isnan(y));
end

%% Function: return image contribution
function dampingContrib = imageContrib(t, L, c, lagDamp, ...
    tVector, xVector)
dampingContrib = 0;

% Right images
imageSign = 1;
tImage = t - 2*L/c;
while tImage >= 0
    nextImage = lagDamp*imageSign*...
        sampleHistory(tVector, xVector(2,:), tImage);
    assert(isreal(nextImage))
    dampingContrib = dampingContrib + nextImage;
    
    % Compute next term params
    tImage = tImage - 2*L/c;
    %         imageSign = -imageSign;
end

% Left images
imageSign = 1;
tImage = t - 2*L/c;
while tImage >= 0
    nextImage = lagDamp*imageSign*...
        sampleHistory(tVector, xVector(2,:), tImage);
    dampingContrib = dampingContrib + nextImage;
    
    % Compute next term params
    tImage = tImage - 2*L/c;
    %         imageSign = -imageSign;
end
end