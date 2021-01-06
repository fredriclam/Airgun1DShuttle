function temperatureRatio = temperatureMachFunction(gamma, M)
% Ratio T/T0 (T0: stagnation temperature) for isentropic flow
temperatureRatio = 1./ (1 + 0.5*(gamma-1)* M .^2 );