gamma = 1.4;

%% Area mach number relation for isentropic flow
% Isentropic flow: analytic relation
isentropic_area = @(M2) ((gamma + 1)/2)^(-(gamma+1)/2/(gamma-1)) ...
    * (1 + 0.5*(gamma-1) * M2).^((gamma+1)/(2*(gamma-1))) ./ sqrt(M2);

%% Differential relation: dT_t/dA, greedy estimate
% Secant estimate for stagnation temperature increase over area change
dTtdA = -0.0035;

polytropic_RHS_inv = @(A, M2, dTtdA) ...
1 ./ ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2) .* ...
(-1 ./ A + 0.5*(1+gamma*M2) * dTtdA);

[A_sub, M2_sub] = ode45(@(A, M2) polytropic_RHS_inv(A, M2, dTtdA), ...
        [1, 12], 1-1e4*eps);
[A_super, M2_super] = ode45(@(A, M2) polytropic_RHS_inv(A, M2, dTtdA), ...
        [1, 8], 1+1e4*eps);
    
% Concatenate sub- and super-sonic, sort according to M^2 values
A = [A_sub; A_super];
M2 = [M2_sub; M2_super];
[M2, perm] = sort(M2);
A = A(perm);

% Plot isentropic flow relation
figure(1); clf;
plot(sqrt(M2), isentropic_area(M2));
hold on
% Plot quadrature result
plot(sqrt(M2), A, '.');

%% Differential relation: dT_t/dA, conservative estimate
% Secant estimate for stagnation temperature increase over area change
dTtdA = -0.035;

polytropic_RHS_inv = @(A, M2, dTtdA) ...
1 ./ ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2) .* ...
(-1 ./ A + 0.5*(1+gamma*M2) * dTtdA);

[A_sub, M2_sub] = ode45(@(A, M2) polytropic_RHS_inv(A, M2, dTtdA), ...
        [1, 12], 1-1e6*eps);
[A_super, M2_super] = ode45(@(A, M2) polytropic_RHS_inv(A, M2, dTtdA), ...
        [1, 8], 1+1e6*eps);
    
% Concatenate sub- and super-sonic, sort according to M^2 values
A = [A_sub; A_super];
M2 = [M2_sub; M2_super];
[M2, perm] = sort(M2);
A = A(perm);

% Append to plot
plot(sqrt(M2), A, '.');

%% Finish plot
xlabel 'M'
ylabel 'A/A^*'
title 'Flow curves'
legend({'Isentropic, analytic', ...
    'Energy loss, quadrature', ...
    'Conservative energy loss, quadrature'})