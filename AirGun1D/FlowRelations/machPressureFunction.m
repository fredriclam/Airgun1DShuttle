function M = machPressureFunction(gamma, pRelative)
if pRelative > 1
    error('Cannot have higher pressure than stagnation.')
end
% From ratio p/p0 (p0: stagnation pressure) for isentropic flow
M = sqrt((pRelative^((1-gamma)/gamma) - 1) * 2 / (gamma - 1));