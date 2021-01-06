gamma = 1.4;

%% Area mach number relation for isentropic flow
% Isentropic flow: analytic relation
isentropic_area = @(M2) ((gamma + 1)/2)^(-(gamma+1)/2/(gamma-1)) ...
    * (1 + 0.5*(gamma-1) * M2).^((gamma+1)/(2*(gamma-1))) ./ sqrt(M2);
% Differential relation
isentropic_RHS = @(M2, A) - A ...
.* ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2);
%% Integrate ODE dA/dM2 = RHS
[M2, A] = ode45(isentropic_RHS, [0.01, 9], 1);
% Scale w.r.t sonic area
A_sonic = pchip(M2, A, 1);
A = A / A_sonic;
%% Plotting isentropic flow relation
figure(1); clf;
plot(sqrt(M2), isentropic_area(M2));
hold on
plot(sqrt(M2),A, '.');
xlabel 'M'
ylabel 'A/A^*'
title 'Isentropic flow'
legend({'Isentropic, analytic', 'Isentropic, quadrature'})

%% Somehow polytropic flow with stagnation temperature change
% Averaged rate of stagnation temperature change per unit M^2 change
dTdM2_vector = [-1, -1e1, -1e2];
% Stagnation temperature at starting state [K]
Tt = 300;

% Differential relation
polytropic_RHS = @(M2, A, dTdM2) - A ...
.* ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2) ...
+ A .* 0.5 *(1 + gamma * M2) * dTdM2 / Tt;

%% Prepare plot for flow relation
figure(2); clf;
plot(sqrt(M2), isentropic_area(M2));
hold on

%% Integrate for each parameter value
for i = 1:length(dTdM2_vector)
    %% Integrate ODE dA/dM2 = RHS
    [M2, A] = ode45(@(M2, A)polytropic_RHS(M2, A, dTdM2_vector(i)), ...
        [0.01, 9], 1);
    % Scale w.r.t sonic area
    A_sonic = pchip(M2, A, 1);
    A = A / A_sonic;
    plot(sqrt(M2),A, '.');
end

%% Finish plotting
xlabel 'M'
ylabel 'A/A^*'
title 'Area change + cooled flow'
legend({'Isentropic, analytic', 'dT_t/dM^2 = -1', 'dT_t/dM^2 = -10', 'dT_t/dM^2 = -100'})