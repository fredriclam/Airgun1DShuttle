% Test data
p = 1e6;
u = 10;
T = 300;
cv = 718;
R = 278;
gamma = 1 + R/cv;
rho = p / (R * T);
e = rho*cv*T + 0.5 * rho * u^2;
sigma = sqrt(2*(gamma-1));

c = @(rho, u, e) sqrt(gamma * R* (e - 0.5 * rho * u^2) / (rho * cv));
Tmat = @(rho, u, e) ...
    [sigma*rho, rho, rho;
     sigma*rho*u, rho*(u+c(rho, u, e)), rho*(u-c(rho, u, e));
     sigma*rho*u.^2/2, ...
       e + (gamma-1)*(e-rho*u.^2/2) + rho*u*c(rho, u, e), ...
       e + (gamma-1)*(e-rho*u.^2/2) - rho*u*c(rho, u, e)]; 
q = @(rho, u, e) [rho; rho*u; e];
char_vec = @(rho, u, e) inv(Tmat(rho, u, e))*q(rho, u, e);

% Entropy-like quantity
s = @(rho, p) p * rho^(-gamma);

% Flux jacobian (verified by eigenvalues)
A = @(rho, u, e) ...
    [ 0, 1, 0;
      0.5*(gamma-3)*u^2, (3-gamma)*u, gamma - 1;
      (gamma-1)*u^3 - gamma*e/rho*u, ...
         gamma * e / rho - 3/2*(gamma-1)*u^2, gamma*u];

T_inst = Tmat(rho, u, e);
A_inst = A(rho, u, e);
[X, L] = eig(A_inst);

prop1 = 2/(gamma-1) * c(rho, u, e) + u;
prop2 = 2/(gamma-1) * c(rho, u, e) - u;
prop3 = s(rho, p);

% Get normalized eigenvector
for i = 1:3
    T_inst_norm(:,i) = T_inst(:,i) ./ norm(T_inst(:,i));
end

% inv(T_inst_norm) * char_vec(rho, u, e) ./ [prop3; prop1; prop2]
(T_inst_norm)' * char_vec(rho, u, e) ./ [prop3; prop1; prop2]