function [dz, p_rear, p_front, p_Mid, subsystemState] = shuttleEvolve( ...
    z, p_L, physConst, chamberSet)
% Computes evolution of the state vector
%   z = [pos; vel; m_L; E_L; m_R; E_R]
% for the shuttle assembly and the operating chamber partitions (denoted
% here L and R). Models momentum on the shuttle, isentropic compression of
% the fluid in the partition and isentropic flow between the partitions.
%
%   Input:
%     z = [pos; vel; m_L; E_L; m_R; E_R]
%     p_L = pressure on left side of shuttle assembly (from firing chamber)
%     physConst = struct of physical constants
%     chamberSet = object of type Chambers
%   Output: dz = [dpos; dvel; ...; dE_R] (time derivatives of state z)
%     p_rear = pressure at rear of shuttle, in operating chamber
%     p_front = pressure at front of shuttle, in operating chamber
%     p_Mid = pressure in middle chamber (if vented to ambient, constant)
%     subsystemState = human-readable struct with all relevant states and
%       differentials

%% Compute thermodynamic state and define helper variables
% Unpack state
m_rear = z(3);
E_rear = z(4);
m_front = z(5);
E_front = z(6);

% Compute elastic penalty for constraint z(1) >= 0
penaltyForce = physConst.shuttleBdryPenaltyStrength * ...
    (z(1) < 0) * abs(z(1));
A_L = physConst.shuttle_area_left;
A_R = physConst.shuttle_area_right;

% Compute mass flow rates
% Geometrically constrained density
rho_rear = m_rear / chamberSet.rearVolume(z(1));
% Energy-constrained temperature
T_rear = E_rear / (m_rear * physConst.c_v);
% Ideal gas pressure
p_rear = rho_rear * physConst.Q * T_rear;
% Geometrically constrained density
rho_front = m_front / chamberSet.frontVolume(z(1));
% Energy-constrained temperature
T_front = E_front / (m_front * physConst.c_v);
% Ideal gas pressure
p_front = rho_front * physConst.Q * T_front;

% In middle chamber: assumptions
midChamberLength = 3.112 * 0.0254;
A_Mid = physConst.shuttle_area_right_rear;

% IMPORTANT: assuming depth of 10m, and also 1000kg/m^3 water
% Atmospheric plus water depth pressure
pAmbient = 1e5 + 9.8*1e3 * 10;
% Initial mid chamber pressure
p0_Mid = pAmbient;

if strcmpi('limit-vented', chamberSet.midChamberMode)
    % p0_Mid = physConst.p_R0;
    % EAGE model: venting valve.
    % Here we assume the equilibration is instantaneous (to be replaced by
    % possible gas exchange to back chamber and exterior)    
    p_Mid = p0_Mid;
elseif strcmpi('limit-closed', chamberSet.midChamberMode)
    % Allowing adiabatic compression
    p_Mid = p0_Mid * (midChamberLength / ...
            (midChamberLength-z(1)))^physConst.gamma;
else
    warning('Mid chamber mode unset')
end

% Define local aliases
g = physConst.gamma;
c_p = physConst.c_v + physConst.Q;

%% Compute flow between partition of operating chamber
% Compute direction of flow and pressure ratio (positive if rear to front)
flowSignL2R = sign(p_rear - p_front);
pRatio = min([p_rear/p_front, p_front/p_rear]);
% Compute mach number due to expansion capping to choked
M = min([1, machPressureFunction(g, pRatio)]);
% Compute stagnation properties
if flowSignL2R >= 0
    pMax = p_rear;
    rhoMax = rho_rear;
    TMax = T_rear;
else
    pMax = p_front;
    rhoMax = rho_front;
    TMax = T_front;
end

flowL2R = ...
    flowSignL2R * chamberSet.gapArea(z(1)) * M * ...
    sqrt(g * pMax * rhoMax) * ...
    (1 + (g-1)/2 * M^2)^(0.5*(-g-1)/(g-1));

%% Compute dz

netForce = ...
      p_L*A_L ...
      - p_Mid * A_Mid ...
      + p_rear*physConst.shuttle_area_right_rear ...
      - p_front*A_R ...
      + penaltyForce;

dz = [z(2);
      netForce/physConst.shuttleAssemblyMass;
      -flowL2R;
      -c_p * flowL2R * TMax - p_rear * physConst.shuttle_area_right_rear * z(2);
      flowL2R;
      c_p * flowL2R * TMax + p_front * A_R * z(2);];

% Create human-readable object with shuttle and chamber data
subsystemState = struct('opChamberRear_m', m_rear, ...
                        'opChamberRear_E', E_rear, ...
                        'opChamberRear_rho', rho_rear, ...
                        'opChamberRear_T', T_rear, ...
                        'opChamberRear_p', p_rear, ...
                        'opChamberFront_m', m_front, ...
                        'opChamberFront_E', E_front, ...
                        'opChamberFront_rho', rho_front, ...
                        'opChamberFront_T', T_front, ...
                        'opChamberFront_p', p_front, ...
                        'penaltyForce', penaltyForce, ...
                        'opChamberFront_area', A_L, ...
                        'opChamberRear_area', A_R, ...
                        'midChamber_length', midChamberLength, ...
                        'midChamber_area', A_Mid, ...
                        'midChamber_p', p_Mid, ...
                        'opChamberFlow_M', M, ...
                        'opChamberFlow_massL2R', flowL2R, ...
                        'opChamberFlow_denergydtL', dz(4), ...
                        'opChamberFlow_denergydtR', dz(6), ...
                        'shuttle_position', z(1), ...
                        'shuttle_velocity', z(2), ...
                        'shuttle_acceleration', dz(2), ...
                        'shuttle_force', netForce, ...
                        'physConst', physConst);
return

%% Legacy damping forces
% The following functions can be used to provide additional sources of
% damping. These sources are typically insignificant compared to the
% chamber gas dynamics.

% Arbitrary damping model
% Input: state vector z --- [pos; vel]
function dampingForce = linearDampingForce(z)
dampingForce = -5e2 * z(2);

% Constant damping force
% Coefficient of friction (would be empirical).
% Input: state vector z --- [pos; vel]
function dampingForce = coulombDampingForce(z, physConst)
magnitude = 0.3 * 9.8 * physConst.shuttleAssemblyMass;
dampingForce = -sign(z(2)) * magnitude;

function dampingForce = constantDampingForce(z)
dampingForce = -10e3; % 2 kN (overpredicting the incompressible turbulent 1/7-law)

% Viscoelastic model: is this just linear damping?
% Input: state vector z --- [pos; vel]
function dampingForce = viscoelasticDampingForce(z)
dampingForce = - 1e3 * z(2);
