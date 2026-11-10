function [dx, apply, P_post, innov, S_vec, sigma_p] = ...
            eskf_update(P_prev, p_ins, v_ins, z_enu, k, gnss_valid, z_vel, eskf_prm)
%ESKF_UPDATE LC error-state KF update with GNSS position + velocity.
%
%     z = [z_pos; z_vel],   y = z - [p_ins; v_ins],   H = [I2 0 0; 0 I2 0]
%
%   P_prev      (8x8)   Covariance predicted to this tick
%   p_ins       (2x1)   INS position, delayed state x_k-  (x_prev(1:2))
%   z_enu       (2x1)   GNSS position fix, ENU [m]
%   k           (1x1)   Degradation multiplier (scales R)
%   gnss_valid  (1x1)   True on exactly one tick per fix
%   v_ins       (2x1)   INS velocity, delayed state x_k-  (x_prev(3:4))
%   z_vel       (2x1)   GNSS velocity fix, ENU [m/s]
%   eskf_prm            R_sigma, Rv_sigma
%
%#codegen

    P = 0.5*(P_prev + P_prev.');

    dx    = zeros(8,1);
    apply = false;
    innov = zeros(2,1);
    S_vec = zeros(3,1);

    if gnss_valid
        H  = [eye(2),     zeros(2,6);
              zeros(2,2), eye(2), zeros(2,4)];                  % 4x8
        Rm = k^2 * blkdiag(eskf_prm.R_sigma^2  * eye(2), ...
                           eskf_prm.Rv_sigma^2 * eye(2));       % 4x4

        y = [reshape(z_enu(1:2), 2, 1) - reshape(p_ins(1:2), 2, 1);
             reshape(z_vel(1:2), 2, 1) - reshape(v_ins(1:2), 2, 1)];   % 4x1
        S = H*P*H.' + Rm;        % 4x4
        K = (P*H.') / S;         % 8x4
        dx = K*y;

        IKH = eye(8) - K*H;
        P   = IKH*P*IKH.' + K*Rm*K.';      % Joseph form
        P   = 0.5*(P + P.');

        apply = true;
        innov = y(1:2);
        S_vec = [S(1,1); S(2,2); S(1,2)];
    end

    P_post  = P;
    sigma_p = sqrt(max(P(1,1) + P(2,2), 0));
end