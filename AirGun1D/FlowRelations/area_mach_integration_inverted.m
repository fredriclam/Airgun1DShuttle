% Area mach number relation

gamma = 1.4;

% Isentropic flow: analytic relation
isentropic_area = @(M2) ((gamma + 1)/2)^(-(gamma+1)/2/(gamma-1)) ...
    * (1 + 0.5*(gamma-1) * M2).^((gamma+1)/(2*(gamma-1))) ./ sqrt(M2);
% isentropic_area_inv = @(A) fsolve(@(M) 

% Integration relation
isentropic_RHS = @(A, M2) -1 ./ A ...
./ ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2);
isentropic_RHS_reverse = @(A_rev, M2) +1 ./ (1-A_rev) ...
./ ((1 - M2) ./ (1 + ((gamma-1)/2) * M2) .* 0.5 ./ M2);

% Integrate ODE dM2/dA = RHS
% Forward:
[t1, y1] = ode45(isentropic_RHS, [1, 3], 1+1e4*eps);
% Backward:
[t2, y2] = ode15s(@(Arev, M2) isentropic_RHS_reverse(Arev,M2), ...
    [0, 0.9], 1-1e4*eps);
% Concat
t = [1-t2; t1];
y = [y2; y1];
% Sort forward
[t, perm] = sort(t);
y = y(perm);

%% Plotting
figure(1); clf;
plot(y,t, '.');
hold on
plot(y, isentropic_area(y));
