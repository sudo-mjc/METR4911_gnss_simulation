function [z_deg, k, gnss_valid_eff, zv_deg] = ...
            gnss_degrade(t, z_pos, p_true, gnss_valid, z_vel, v_true, deg_prm)
%GNSS_DEGRADE  Time-scheduled GNSS degradation, position AND velocity.
%
%     z_deg  = p_true + k*(z_pos - p_true)     position error std grows by k
%     zv_deg = v_true + k*(z_vel - v_true)     velocity error std grows by k
%     k      -> ESKF: BOTH position and velocity R are scaled by k^2
%     gnss_valid_eff = gnss_valid AND (k < k_out)   outage drops both
%
%   Why the same k for velocity: the degraded environments modelled here
%   (canopy / urban canyon) degrade the whole fix - fewer satellites,
%   worse geometry (DOP) and multipath affect Doppler velocity as well as
%   pseudorange position. Doppler usually degrades somewhat less, so a
%   shared k is a slightly conservative (pessimistic) choice. It also needs
%   NO new deg_prm fields, so build_sim_params.m is unchanged for this block.
%   Because R is scaled by the same k, the filter stays consistent with the
%   error actually injected.
%
%   Inputs
%     t          : simulation time (Clock block)                     [s]
%     z_pos      : GNSS position in the local ENU frame              [m]
%     p_true     : truth position (3D padding fine)                  [m]
%     gnss_valid : fix-available pulse from gnss_schedule
%     z_vel      : GPS block Velocity output, local ENU (1x3)        [m/s]
%                  (GPS block 'Reference frame' must be ENU - if NED,
%                   the first two elements are [vN vE], not [vE vN])
%     v_true     : truth velocity, same frame (truth.velocity)       [m/s]
%     deg_prm    : flat struct from build_sim_params (Parameter scope)
%                  en (0/1), t (1xN segment starts), k (1xN), k_out
%
%   Outputs
%     z_deg, k, gnss_valid_eff
%     zv_deg     : (2x1) degraded horizontal GNSS velocity [vE; vN]  [m/s]
%
%   With deg_prm.en = 0: z_deg = z_pos, zv_deg = z_vel(1:2), k = 1.
%#codegen

% Position: truth reshaped to the size of z_pos (fixed size on every path)
p = reshape(p_true(1:numel(z_pos)), size(z_pos));

% Velocity: horizontal channel only, always a 2x1 column for the ESKF
zv = reshape(z_vel(1:2),  2, 1);
vt = reshape(v_true(1:2), 2, 1);

k = 1;
if deg_prm.en ~= 0
    % Piecewise-constant lookup: k of the last segment that has started
    for i = 1:numel(deg_prm.t)
        if t >= deg_prm.t(i)
            k = deg_prm.k(i);
        end
    end
end

z_deg  = p  + k*(z_pos - p);       % k = 1 -> exact pass-through
zv_deg = vt + k*(zv - vt);
gnss_valid_eff = logical(gnss_valid) && (k < deg_prm.k_out);
end