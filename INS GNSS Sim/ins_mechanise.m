function [x_next, x_nav, f_b, aNav, w_corr] = ins_mechanise(x_prev, dx, apply, ...
                                                    accelRead, gyroRead, ins_prm)
%INS_MECHANISE Correct-then-propagate strapdown mechanisation (2D, ENU).
%
%   x = [pE pN vE vN psi b_ax b_ay b_gz]'
%
%   x_prev     (8x1)  state from the Unit Delay
%   dx         (8x1)  correction from the ESKF
%   apply      (1x1)  true on update ticks
%   accelRead  (1x3)  raw accelerometer, body [m/s^2]
%   gyroRead   (1x3)  raw gyro, body [rad/s]
%   ins_prm               dt, accelInv, feedback
%
%   x_next     (8x1)  -> Unit Delay
%   x_nav      (8x1)  navigation output at THIS instant (log this)
%   f_b        (2x1)  bias-corrected specific force -> ESKF
%   aNav       (2x1)  nav-frame acceleration
%#codegen

    x = x_prev;

    %% 1. Correct
    if apply && ins_prm.feedback
        x    = x + dx;
        x(5) = wrap2pi(x(5));
    end
    x_nav = x;                       % output BEFORE propagating: refers to t_k

    %% 2. Propagate
    p   = x(1:2);
    v   = x(3:4);
    psi = x(5);
    b_a = x(6:7);
    b_g = x(8);

    a_meas = reshape(accelRead(1:2), 2, 1);       % force column, whatever the bus gives
    f_b    = ins_prm.accelInv * a_meas - b_a;

    c = cos(psi);
    s = sin(psi);
    aNav = [c, -s; s, c] * f_b;                   % body -> nav

    p_n   = p + v*ins_prm.dt + 0.5*aNav*ins_prm.dt^2;
    v_n   = v + aNav*ins_prm.dt;
    psi_n = wrap2pi(psi + (gyroRead(3) - b_g)*ins_prm.dt);

    w_corr = gyroRead(3) - b_g;      % 1x1

    x_next = [p_n; v_n; psi_n; b_a; b_g];
end