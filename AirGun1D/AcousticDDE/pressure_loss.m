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

cushLength = 0.5 * 0.0254;      % [m]
accelLength = 2.5 * 0.0254;     % [m]
pRFactor = @(xi) (xi <= accelLength) * (1) + ...
    (xi > accelLength) * ...
    (cushLength / (cushLength - (xi - accelLength))) .^ gamma; % [-]
pR = @(xi) p0 * pRFactor(xi); % [Pa]

A_L = 10* 16 * (0.0254^2);        % [m^2] -- from test_launch_script test case 1
A_R = 0.0506/0.0586 * A_L;    % [m^2]

% L = c * pi / sqrt((A_R * p0 * gamma / cushLength)/m); % Resonant length [m]
% 
% f0 = 190; % [Hz] from FFT, short chamber
% L = c / (2*f0);

L = 0.6; % physical measurement of short airgun

%% Compute equilibrium position
forceFactor = ((A_L * p0) / (A_R * p0))^(1/gamma);
xi_eq = (forceFactor-1)/forceFactor * cushLength;

%% Solve and plot
figure(1); clf;
tFinal = 0.025;

% Solve
[tVector1, xVector1] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.01, 0);
subplot(1,3,1);
plot(tVector1, xVector1(1,:));
hold on
[tVector1, xVector1] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.01, 1e8);
plot(tVector1, xVector1(1,:));
title('L = 0.6 m; depressurization timescale 0.01s')

% Solve
[tVector2, xVector2] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.03, 0);
subplot(1,3,2);
plot(tVector2, xVector2(1,:));
hold on
[tVector2, xVector2] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.03, 1e8);
plot(tVector2, xVector2(1,:));
title('L = 0.6 m; depressurization timescale 0.03s')

% Solve
[tVector3, xVector3] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.1, 0);
subplot(1,3,3);
plot(tVector3, xVector3(1,:));
hold on
[tVector3, xVector3] = DDESolve(tFinal, true, true, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, ...
    A_R, A_L, 0.1, 1e8);
plot(tVector3, xVector3(1,:));
title('L = 0.6 m; depressurization timescale 0.1s')

for i = 1:3
    subplot(1,3,i);
    xlabel('time [s]')
    ylabel('shuttle position [m]')
    ylim((accelLength + cushLength) *[-.1 1]);
    hold on
    plot([0, tFinal], accelLength*[1 1], '--');
end

% hold on;
% subplot(2,3,2);
% plot(tVector1, xVector1(2,:));
% hold on;
% title('Velocity [m/s] vs. time [s]')

% subplot(2,3,4);
% dt = tVector1(2)-tVector1(1);
% [f1, X1] = outil.spect(xVector1(1,:), dt);
% loglog(f1, X1);

%% Function: DDE solve with fixed timestep
function [tVector, xVector] = DDESolve(tFinal, globalDamp, lagDamp, ...
    m, L, p0, T0, gamma, R, c, rho0, cushLength, accelLength, pR, A_R, ...
    A_L, depressurizationTimescale, penaltySpringK)

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
    
    % Artificial decay term
    pamb = 1e5;
    p_L = pamb + (p0-pamb) * exp(-t/depressurizationTimescale);
    
    % Sum contributions over domain of influence:
    damping = x(2) + ...
        imageContrib(t, L, c, lagDamp, ...
        tVector, xVector);
    
    % Build RHS at time t
    dx = [x(2); ...
          1/m * (...
              A_L * p_L - A_R * pR(x(1)) ...
              - globalDamp * rho0*c * A_L * damping) ...
              - penaltySpringK * (x(1) < 0) * x(1) 
         ] * dt;
    
    xGuess = x + dx / 2;
    
    % Build approximate RHS at time t + dt/2; use same xVector (stability
    % condition needs to be respected)
    damping = xGuess(2) + ...
        imageContrib(t + dt/2, L, c, lagDamp, ...
        tVector, xVector);
    
    dx = [xGuess(2); ...
          1/m * (...
              A_L * p_L - A_R * pR(xGuess(1)) ...
              - globalDamp * rho0*c * A_L * damping) ...
              - penaltySpringK * (x(1) < 0) * x(1)
         ] * dt;
    % Get next x from Forward Euler operator
    xVector(:, i+1) = x + dx;
    
%     % Hard enforement
%     if xVector(1, i+1) < 0
%         xVector(1, i+1) = 0;
%         if xVector(2, i+1) < 0
%             xVector(2, i+1) = -xVector(2, i+1);
%         end
%     end
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