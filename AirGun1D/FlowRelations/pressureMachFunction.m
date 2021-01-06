function pressureRatio = pressureMachFunction(gamma, M)
% Ratio p/p0 (p0: stagnation pressure) for isentropic flow
pressureRatio = (1 + 0.5*(gamma-1)* M .^2 ).^ ...
                (-gamma/(gamma-1));