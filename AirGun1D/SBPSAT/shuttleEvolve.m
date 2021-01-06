function [dz, p_rear, p_front, p_Mid] = shuttleEvolve(z, p_L, physConst, opChamber)
% Computes evolution of the state vector
%   z = [pos; vel; m_L; E_L; m_R; E_R]
% for the shuttle assembly and the operating chamber partitions (denoted
% here L and R). Models momentum on the shuttle, isentropic compression of
% the fluid in the partition and isentropic flow between the partitions.
%
%   Input:
%     z = [pos; vel; m_L; E_L; m_R; E_R]
%     p_L = pressure on left side of shuttle assembly (from firing chamber)
%     physConst = struct of physical constnats
%     opChamber = object of type OperatingChamber
%   Output: dz = [dpos; dvel; ...; dE_R] (time derivatives of state z)
% 
% Last update: 2020-09-23 (minor fix and cleanup).

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
rho_rear = m_rear / opChamber.rearVolume(z(1));
% Energy-constrained temperature
T_rear = E_rear / (m_rear * physConst.c_v);
% Ideal gas pressure
p_rear = rho_rear * physConst.Q * T_rear;
% Geometrically constrained density
rho_front = m_front / opChamber.frontVolume(z(1));
% Energy-constrained temperature
T_front = E_front / (m_front * physConst.c_v);
% Ideal gas pressure
p_front = rho_front * physConst.Q * T_front;

% In middle chamber: assumptions
midChamberLength = 3.112 * 0.0254;
A_Mid = physConst.shuttle_area_right_rear;
p0_Mid = physConst.p_R0;
% EAGE model: venting valve.
% Here we assume the equilibration is instantaneous (to be replaced by
% possible gas exchange to back chamber and exterior)
p0_Mid = 9.8*1e4; % IMPORTANT: assuming depth of 10m, and also 1000kg/m^3
% Allowing adiabatic compression
% p_Mid = p0_Mid * (midChamberLength / ...
%         (midChamberLength-z(1)))^physConst.gamma;
% Alternative assumption: instantaneous equilibration of mid chamber
% pressure
p_Mid = p0_Mid;

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

efficiencyFactor = 1.0;

flowL2R = efficiencyFactor * ...
    flowSignL2R * opChamber.gapArea(z(1)) * M * ...
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

% assert(all(isreal(dz)) && all(~isnan(dz)));

return

%% Legacy damping forces
% The following functions can be used to provide additional sources of
% damping.

% Arbitrary damping model
% Input: state vector z --- [pos; vel]
function dampingForce = linearDampingForce(z)
dampingForce = -5e2 * z(2);

% Constant damping force
% Coefficient of friction (would be empirical). Doesn't really do much!
% Input: state vector z --- [pos; vel]
function dampingForce = coulombDampingForce(z, physConst)
magnitude = 0.3 * 9.8 * physConst.shuttleAssemblyMass;
dampingForce = -sign(z(2)) * magnitude;

function dampingForce = constantDampingForce(z)
dampingForce = -10e3; % 2 kN (overpredicting the INcompressible turbulent 1/7-law)

% Viscoelastic model: is this just linear damping?
% Input: state vector z --- [pos; vel]
function dampingForce = viscoelasticDampingForce(z)
dampingForce = - 1e3 * z(2);
