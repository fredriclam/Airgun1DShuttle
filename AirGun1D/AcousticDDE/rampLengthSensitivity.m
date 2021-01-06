cushLength = 0.5 * 0.0254;      % [m]
accelLength = 2.5 * 0.0254;     % [m]

rampLengths = linspace(0.001, 1, 1000) * accelLength;
critTaus = nan(size(rampLengths));

for i = 1:length(rampLengths)
    rampLength = rampLengths(i);
    critTaus(i) = critical_depressurization_timescale(rampLength);
end

figure(1); clf;
subplot(1,2,1);
plot(rampLengths, critTaus);
xlabel('\eta_{ramp}')
ylabel('\tau_{crit}')

subplot(1,2,2);
loglog(rampLengths, critTaus);
xlabel('\eta_{ramp}')
ylabel('\tau_{crit}')

%% Critical tau is near where c_2 == 0 (for nu2, exponential growth mode)
function tau_star = critical_depressurization_timescale(rampLength)
% Fix parameters characteristic to airgun
p0 = 13609200;
pR0 = p0;
pL0 = p0;
c = 347.1887;
rho0 = 158.0627;
A_Rbehind = 0.0802;
A_R = 0.0891;
A_L = 0.1032;
% rampLength = 0.0032;
m = 28.6020;

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

% Approximate the critical tau with the exponential 
tau_star = fzero(@(tau) [0 1] * mode_weights(tau), [0, 1e-2]);
end