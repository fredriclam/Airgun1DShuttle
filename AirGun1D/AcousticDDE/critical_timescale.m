% Saved parameters
p0 = 13609200;
pR0 = p0;
pL0 = p0;
c = 347.1887;
rho0 = 158.0627;
A_Rbehind = 0.0802;
A_R = 0.0891;
A_L = 0.1032;
rampLength = 0.0032;
m = 28.6020;

%% First attempt for critical tau
% sigma = 1/sqrt(pR0 * A_Rbehind / rampLength /m);
% alpha = @ (tau) -tau^2 * pL0 * A_L / ...
%     (pR0 * A_Rbehind / rampLength * tau^2 + rho0 * c * A_L * tau - m);
% 
% group = @(tau) -(pR0 * A_Rbehind/rampLength * tau^2 + rho0 * c * A_L * tau - m) / (tau^2 * pL0 * A_L);
% f = @(tau) sqrt(pR0*A_Rbehind / rampLength /m) - tau * (1 + A_R/A_Rbehind * rampLength * group(tau));
% tau_star = fzero(f, 1e-4)

%% Parametric
a1 = m;
a2 = rho0 * c * A_L;
a3 = pR0 * A_Rbehind / rampLength;
b1 = pR0 * A_R;
b2 = pL0 * A_L;

% Imag frequencies of system response
nu1 = (-a2 - sqrt(a2^2 + 4*a1*a3)) / (2*a1);
nu2 = (-a2 + sqrt(a2^2 + 4*a1*a3)) / (2*a1);

% System response constants
BC_particular = @(tau) [b1/a3 + b2*tau^2 / (a1 - a2*tau - a3*tau^2);
                        -b2*tau / (a1 - a2*tau - a3*tau^2)];
mode_weights = @(tau) -[1, 1; nu1, nu2] \ BC_particular(tau);
% Modes
response_vector = @(t) [exp(nu1*t); ...
                        exp(nu2*t)];

% For some tau, plot the asymptotic behaviour
tau = 1.8831e-3; % Critical tau
tau = 1.8831e-3;
response = @(t) mode_weights(tau)' * response_vector(t) + b1/a3 + ...
    b2*tau^2 / (a1 - a2*tau - a3*tau^2) * exp(-t/tau);
ezplot(response, [0, 5e-3])
% ylim(rampLength*[-1, 1]);

%% Critical tau is near where c_2 (for nu2, exponential growth mode)
figure(1); clf;
tFinal = 0.015;
tau_star = fzero(@(tau) [0 1] * mode_weights(tau), 2e-3);
critical_weights = mode_weights(tau_star);
response = @(t) [critical_weights(1) 0] * response_vector(t) + b1/a3 + ...
b2*tau^2 / (a1 - a2*tau - a3*tau^2) * exp(-t/tau);
ezplot(response, [0, tFinal])
hold on

% Negative eps exponential
response = @(t) [critical_weights(1) -1e5*eps] * response_vector(t) + b1/a3 + ...
b2*tau^2 / (a1 - a2*tau - a3*tau^2) * exp(-t/tau);
ezplot(response, [0, tFinal])

% Positive eps exponential
response = @(t) [critical_weights(1) eps] * response_vector(t) + b1/a3 + ...
b2*tau^2 / (a1 - a2*tau - a3*tau^2) * exp(-t/tau);
ezplot(response, [0, tFinal])

ylim([0, 5e-3])
plot([0, 0.05], rampLength*[1 1])
legend({'Critical tau: zero positive exponential', ...
    'negative eps of exponential',...
    'positive eps of exponential','Ramp length'}, 'location', 'best')