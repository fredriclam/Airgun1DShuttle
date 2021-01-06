% Text file with code work in progress

% Compute sonic area as a function of conserved variables
q = [_, _, _];
% Temperature


p_R = schm.p(q_R);
rho_R = q_R(1);
u_R =  q_R(2)/q_R(1);
T_R = p_R / rho_R / physConst.Q;
c_R = sqrt(physConst.gamma * physConst.Q * T_R);
M_R = u_R / c_R;

% Private alias
gamma_ = physConst.gamma;
% Compute ratio of area to sonic area
A_on_Asonic_ratio = ...
    ((gamma_+1)/2)^(-(gamma_+1)/2/(gamma_-1)) * ...
    (1 + (gamma_-1)/2 * M_R^2 )^((gamma_+1)/2/(gamma_-1)) ./ M_R;
% Compute sonic area
A_sonic = A_cs / A_on_Asonic_ratio;

% Isentropic ratios w.r.t. reference state
area_ratio = @(M) ((gamma_+1)/2)^(-(gamma_+1)/2/(gamma_-1)) * ...
    (1 + (gamma_-1)/2 * M.^2 ).^((gamma_+1)/2/(gamma_-1)) ./ M;
pressure_ratio = @(M) (1 + 0.5*(gamma-1)* M .^2 ) .^ (-gamma/(gamma-1));
temperature_ratio = @(M) 1./ (1 + 0.5*(gamma-1)* M .^2 );

% Inverse problem for finding subsonic mach number at port
obj_func = @(M_port) A_Port / A_cs - sonic_area_ratio_port(M_port) / sonic_area_ratio_R;
M_port = fsolve(obj_func, 0.5);
p_port = pressure_ratio(M_port) / pressure_ratio(M_R) * p_R;
T_port = temperature_ratio(M_port) / temperature_ratio(M_R) * T_R;
rho_port = p_port / physConst.Q / T_port;