function P_next = eskf_predict(P_post, psi_nav, f_b, eskf_prm)
%ESKF_PREDICT Covariance prediction for the error-state KF.
%
%   P_post   (8x8)  covariance after this tick's update
%   psi_nav  (1x1)  CORRECTED INS heading at this tick [rad]
%   f_b      (2x1)  bias-corrected specific force, body frame [m/s^2]
%   eskf_prm        dt, Na, Ng, Ka, Kg, Q_scale
%
%   P_next   (8x8)  covariance predicted to the next tick -> Unit Delay
%#codegen

    dt = eskf_prm.dt;
    c  = cos(psi_nav);
    s  = sin(psi_nav);

    Rbn   = [ c, -s;  s,  c];      % body -> nav,  R(psi)
    Rbn90 = [-s, -c;  c, -s];      % dR/dpsi     = R(psi + pi/2)

    A = zeros(8,8);
    A(1:2,3:4) = eye(2);           % d(dp)/dt = dv
    A(3:4,5)   = Rbn90 * f_b;      % heading error -> velocity error
    A(3:4,6:7) = -Rbn;             % accel bias    -> velocity error
    A(5,8)     = -1;               % gyro bias     -> heading error

    Adt = A*dt;
    F   = eye(8) + Adt + 0.5*(Adt*Adt);

    G = zeros(8,6);
    G(3:4,1:2) = -Rbn;             % accel white noise
    G(5,3)     = -1;               % gyro white noise
    G(6:7,4:5) = eye(2);           % accel bias random walk
    G(8,6)     = 1;                % gyro bias random walk

    Qc = diag([eskf_prm.Na^2, eskf_prm.Na^2, eskf_prm.Ng^2, ...
               eskf_prm.Ka^2, eskf_prm.Ka^2, eskf_prm.Kg^2]) * eskf_prm.Q_scale;

    M  = G*Qc*G.';
    Qd = 0.5*(F*M*F.' + M) * dt;   % trapezoidal (van Loan) discretisation

    P_next = F*P_post*F.' + Qd;
    P_next = 0.5*(P_next + P_next.');
end