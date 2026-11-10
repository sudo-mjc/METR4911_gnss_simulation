function truth = gen_truth_trajectory(cfg)
%GEN_TRUTH_TRAJECTORY Planar truth trajectory with random stops.
%
% How: the original motion formulas are evaluated on a "motion clock"
% tau instead of t. A gate s(t) in [0,1] is 1 while walking, ramps to 0
% for a stop and back to 1 afterwards, and tau = integral of s dt. So
% while stopped tau (and therefore the heading) is frozen.
%   psi(t) = psi_orig(tau)          r(t) = s * r_orig(tau)
%   V(t)   = s * V_orig(tau)        Vdot = sdot*V_orig(tau) + s^2*Vdot_orig(tau)
%
% Stop times come from a local RandStream seeded with cfg.sim.seed:
% repeatable, and the global rng is untouched.
%
% Outputs: identical names, sizes and types to gen_truth_trajectory.m.
%   truth.t, .position, .velocity, .acceleration, .accelBody,
%   .angularVelocity, .psi, .orientation

    dt = cfg.sim.dt_imu;
    N  = round(cfg.sim.T / dt) + 1;
    t  = (0:N-1)' * dt;

%% Stop schedule (random but repeatable)
% ---------------------------------------------------------------
    T_walk = [60 180];   % (s) walking time between stops
    T_stop = [10 30];    % (s) time fully stationary
    T_ramp = 3;          % (s) slow-down / speed-up time

    rs = RandStream('mt19937ar', 'Seed', cfg.sim.seed);
    U  = @(lim) lim(1) + (lim(2) - lim(1)) * rand(rs);

    % Stop i: decelerate from a(i), stationary from a(i)+T_ramp to b(i),
    % accelerate from b(i) to b(i)+T_ramp.
    a = [];  b = [];
    tc = 0;
    while true
        ai = tc + U(T_walk);
        if ai >= t(end), break; end
        bi = ai + T_ramp + U(T_stop);
        a(end+1) = ai;   b(end+1) = bi;   %#ok<AGROW>
        tc = bi + T_ramp;
    end

%% Gate s(t), its derivative, and the motion clock tau(t)
% ---------------------------------------------------------------
% Raised-cosine ramps, so s and sdot are continuous. Each stop removes
% exactly (b - a) seconds of motion time.
    s    = ones(N,1);
    sdot = zeros(N,1);
    lost = zeros(N,1);           % integral of (1 - s)
    w    = pi / T_ramp;

    for i = 1:numel(a)
        u = t - a(i);
        dec = (u >= 0) & (u < T_ramp);                    % slowing down
        s(dec)    = 0.5*(1 + cos(w*u(dec)));
        sdot(dec) = -0.5*w*sin(w*u(dec));
        lost(dec) = lost(dec) + 0.5*(u(dec) - sin(w*u(dec))/w);

        stp = (t >= a(i) + T_ramp) & (t < b(i));         % stationary
        s(stp)    = 0;
        lost(stp) = lost(stp) + (t(stp) - a(i)) - 0.5*T_ramp;

        u = t - b(i);
        acc = (u >= 0) & (u < T_ramp);                    % speeding up
        s(acc)    = 0.5*(1 - cos(w*u(acc)));
        sdot(acc) = 0.5*w*sin(w*u(acc));
        lost(acc) = lost(acc) + (b(i) - a(i)) - 0.5*T_ramp ...
                              + 0.5*(u(acc) + sin(w*u(acc))/w);

        done = u >= T_ramp;                               % stop complete
        lost(done) = lost(done) + (b(i) - a(i));
    end
    tau = t - lost;

%% Motion profile (original formulas, evaluated at tau)
% ---------------------------------------------------------------
    wV  = 2*pi / cfg.traj.T_V;
    Vo  = cfg.traj.V0 + cfg.traj.dV * sin(wV*tau);
    Vdo = cfg.traj.dV * wV .* cos(wV*tau);
    assert(min(Vo) > 0, 'V0/dV combination allows speed to reach zero or go negative');

    A1 = cfg.traj.yawRateAmp(1);   w1 = 2*pi / cfg.traj.yawRatePer(1);
    A2 = cfg.traj.yawRateAmp(2);   w2 = 2*pi / cfg.traj.yawRatePer(2);
    ro  = A1 * sin(w1*tau) + A2 * cos(w2*tau);
    psi = cfg.traj.psi0 + (A1/w1) * (1 - cos(w1*tau)) + (A2/w2) * sin(w2*tau);

    % Apply the gate
    V    = s .* Vo;
    Vdot = sdot .* Vo + s.^2 .* Vdo;
    r    = s .* ro;

%% Generate trajectory data
% ---------------------------------------------------------------
    vE = V .* cos(psi);
    vN = V .* sin(psi);

    aE = Vdot .* cos(psi) - V .* r .* sin(psi);
    aN = Vdot .* sin(psi) + V .* r .* cos(psi);

    aF = Vdot;
    aL = V .* r;

    pE = cfg.traj.p0(1) + cumtrapz(t, vE);
    pN = cfg.traj.p0(2) + cumtrapz(t, vN);

    Z = zeros(N,1);

    truth.t                = t;
    truth.position         = [pE, pN, Z];
    truth.velocity         = [vE, vN, Z];
    truth.acceleration     = [aE, aN, Z];
    truth.accelBody        = [aF, aL, Z];   % F, L, U
    truth.angularVelocity  = [Z, Z, r];     % yaw rate about Up

    truth.psi         = psi;
    truth.orientation = quaternion([cos(psi/2), Z, Z, sin(psi/2)]);
end