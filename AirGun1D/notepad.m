%%
f = precomputeMachAreaFunction(1.4)
f(pi)
%%
A = pi;
g = (0.5*(gamma+1))^(0.5*(-gamma+1)/(gamma-1))

p = 1/pi;
gamma = 1.4;
M = sqrt((p^((1-gamma)/gamma) - 1) * 2 / (gamma - 1));