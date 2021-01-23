qStar0 = fsolve(...
@(q) diag([1e5 1 1]) * [
essentialConstraint(q);
outgoingCharBCConstraints(q);
], q_R, optimoptions('fsolve','FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13));

qStar1 = fsolve(...
    @(q) diag([1 1 1]) * [
essentialConstraint(q);
outgoingCharBCConstraints(q);
], q_R, optimoptions('fsolve','FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13));

qStar2 = fsolve(...
@(q) diag([1 1e5 1e5]) * [
essentialConstraint(q);
outgoingCharBCConstraints(q);
], q_R, optimoptions('fsolve','FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13));