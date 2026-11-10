function cfg = sim_config()
%SIM_CONFIG  Single source of truth for the 2D INS/GNSS-RTK sensor-fusion-toolbox sim.
%
% Every numeric field below maps either directly onto a Simulink block
% dialog parameter or onto the truth-trajectory
% generator.
%

cfg = struct();

%% ---------------------------------------------------------------------
%  Simulation timing
%  ---------------------------------------------------------------------
cfg.sim.dt_imu    = 0.01;       % IMU / model base sample time (s)
cfg.sim.T         = 3600;         % total run length (s)
cfg.sim.seed      = 44;         % Sim RNG seed

%% ---------------------------------------------------------------------
%  Truth trajectory (planar unicycle / coordinated-turn model)
%  ---------------------------------------------------------------------
% Continuously-curving path parameters
cfg.traj.V0        = 0.9;       % nominal walking speed (m/s)
cfg.traj.dV        = 0.3;       % speed oscillation amplitude (m/s)
cfg.traj.T_V        = 37;       % speed oscillation period (s)
cfg.traj.yawRateAmp = [0.045, 0.025];  % yaw-rate harmonic amplitudes (rad/s)
cfg.traj.yawRatePer = [48, 23];        % yaw-rate harmonic periods (s), incommensurate
cfg.traj.psi0       = 0;               % initial heading (rad)
cfg.traj.p0         = [0; 0];          % initial position, ENU north/east (m)

%% ---------------------------------------------------------------------
%  IMU block (Sensor Fusion and Tracking Toolbox / Navigation Toolbox)
%  Values are the ICM-20948 datasheet noise densities.
%  ---------------------------------------------------------------------
cfg.imu.ReferenceFrame  = "ENU";
cfg.imu.Seed            = cfg.sim.seed;
cfg.imu.accelInv        = -1;
%% IMU — ICM-20948 bench characterisation (static, 171 min, 10 Nov 2026)
% Allan-variance measured values (effective, what the filter should see):
%   accel N = 1.614e-3 m/s^2/rtHz   K <= 6.30e-5 m/s^2*rtHz   (worst of x/y)
%   gyro  N = 1.230e-4 rad/s/rtHz   K <= 3.96e-6 rad/s*rtHz   (z)
% Dialog values below are compensated for the IMU block's NoiseType scaling
% (white x sqrt(2), random walk / sqrt(2)) so the block reproduces the above.
cfg.imu.meas.Na = 1.614e-3;   % accel white      (m/s^2)/rtHz  (worst of x/y)
cfg.imu.meas.Ka = 6.30e-5;    % accel bias RW    (m/s^2)/rts   (upper bound)
cfg.imu.meas.Ng = 1.230e-4;   % gyro  white      (rad/s)/rtHz  (z)
cfg.imu.meas.Kg = 3.96e-6;    % gyro  bias RW    (rad/s)/rts   (upper bound)

cfg.imu.blk.white = sqrt(2);
cfg.imu.blk.rw    = 1/sqrt(2);

% Block dialog values
cfg.imu.accel.VelocityRandomWalk     = cfg.imu.meas.Na * cfg.imu.blk.white;
cfg.imu.accel.AccelerationRandomWalk = cfg.imu.meas.Ka * cfg.imu.blk.rw;
cfg.imu.accel.BiasInstability        = 0;
cfg.imu.gyro.AngleRandomWalk         = cfg.imu.meas.Ng * cfg.imu.blk.white;
cfg.imu.gyro.RateRandomWalk          = cfg.imu.meas.Kg * cfg.imu.blk.rw;
cfg.imu.gyro.BiasInstability         = 0;

cfg.imu.accel.ConstantBias = [0.020, -0.015, 0];
cfg.imu.gyro.ConstantBias = [0, 0, deg2rad(0.10)];

% Deliberate initial INS error, as if initialised from one fix + rough compass
cfg.ins.x0_err = [0.2; -0.2; ...       % position: one fix, ~R_sigma (m)
                  0.02; -0.02; ...     % velocity: GNSS velocity at start-up (m/s)
                  deg2rad(5); ...      % heading: compass-like (rad)
                  0; 0; 0];            % biases: unknown, start at zero
cfg.ins.x0 = [0; 0; cfg.traj.V0; 0; cfg.traj.psi0; 0; 0; 0] + ...
                cfg.ins.x0_err;

%% ---------------------------------------------------------------------
%  GPS block (Navigation Toolbox / Sensor Fusion and Tracking Toolbox)
%  ---------------------------------------------------------------------
cfg.gnss.ReferenceFrame        = "ENU";
cfg.gnss.PositionInputFormat   = "Local";
cfg.gnss.ReferenceLocation     = [0, 0, 0];   % [lat lon alt], arbitrary origin -- only the LOCAL
                                              % Cartesian result is used downstream, so the absolute
                                              % geodetic origin is inconsequential for this study.
cfg.gnss.sigma_h = 0.05;                                          % (m) intended per-axis noise (RTK, 5 cm)
cfg.gnss.HorizontalPositionAccuracy = cfg.gnss.sigma_h / sqrt(2); % block outputs sqrt(2)x this per axis (measured)
cfg.gnss.VerticalPositionAccuracy   = 1.2;    % (m)
cfg.gnss.VelocityAccuracy           = 0.2;    % (m/s)
cfg.gnss.DecayFactor                = 0.95;      % 0 for the consistency test, 0.95 for realistic runs
cfg.gnss.Seed                       =cfg.sim.seed;

% --- Duty cycle (NOT modelled by the GPS block)
cfg.gnss.f_gnss   = 10;   % fix rate within an active window (Hz)
cfg.gnss.T_active = 1.0;  % active-window length (s)
cfg.gnss.T_sleep  = 20;  % fixed baseline sleep length (s) (

cfg.gnss.prm.dt = cfg.sim.dt_imu;   % GNSS time step (must match IMU for sim)
cfg.gnss.prm.T_active = 1.0;        % GNSS active period
cfg.gnss.prm.T_sleep = 4.0;         % GNSS sleep period
cfg.gnss.prm.f_gnss = 10;           % GNSS frequency

%% GNSS Signal Degradation
cfg.gnss.deg.en = true;            % Enable flag
%                      open  reduce  open  canyon  open  denial   open
cfg.gnss.deg.t     = [  0    400     700   1200   1500    1800     2600 ];  % (s) segment starts
cfg.gnss.deg.k     = [  1      3     1      6      1      11        1 ];  % (-) error multiplier
cfg.gnss.deg.k_out = 10;            % (-) k at or above this = no usable fix

%% ESKF
cfg.eskf.R_scale  = 1.41;           % R measurement noise scale
cfg.eskf.Q_scale  = 1.0;            % Estimation scale
cfg.eskf.R_sigma  = cfg.gnss.sigma_h;      % m per axis
cfg.eskf.Rv_sigma = cfg.gnss.VelocityAccuracy;
cfg.eskf.feedback = true;       % false reproduces the feed-forward failure mode

% 1-sigma initial uncertainties
cfg.eskf.P0_diag = [ 0.3, 0.3, ...                 % position   (m)
                     0.05, 0.05, ...             % velocity   (m/s)
                     deg2rad(8), ...           % heading    (rad)
                     0.05, 0.05, ...           % accel bias (m/s^2)
                     deg2rad(0.3) ].^2;        % gyro bias  (rad/s)

%% ---------------------------------------------------------------------
%  Adaptive duty-cycle controller
%  ---------------------------------------------------------------------
cfg.ctrl.T_run      = cfg.sim.T;% Run time for testing control functions
cfg.ctrl.T_warmup   = 15;      % (s) receiver held on at start so the ESKF converges
cfg.ctrl.t_reacq    = 4;       % (s) wake -> first fix
cfg.ctrl.R_req      = 1.0;     % (m)
cfg.ctrl.p_req      = 0.85;    % Statistical requirement
cfg.ctrl.sigma_thr  = cfg.ctrl.R_req / sqrt(-log(1 - cfg.ctrl.p_req)); % (m)
cfg.ctrl.sigma_set  = 0.90 * cfg.ctrl.sigma_thr;  % (m) 2.17 m target the controller aims at
cfg.ctrl.open_loop  = true;    % true = hold T_fix, false = closed loop
cfg.ctrl.T_fix      = 10;      % (s) open-loop sleep
cfg.ctrl.T_sleep0   = 10;      % (s) first closed-loop sleep
cfg.ctrl.T_min      = 5;       % (s) sleep limits
cfg.ctrl.T_max      = 60;      % Suggest 35.3s
cfg.ctrl.n_skip     = 5;       % ignore the first coast (it follows the warm-up)
cfg.ctrl.Ki         = 0.5;     % Integral gain on the crossing error
cfg.ctrl.Kp         = 0.3;     % Proportional gain on crossing error

cfg.ctrl.adapt_en   = false;   % Enable boolean for adaptive or bang-bang control
cfg.ctrl.m_fix      = 2.5;     % (-) growth exponent: between velocity (2) and bias (3-4) dominated
cfg.ctrl.dT_default = 1.0;     % (s) Fixed time step when prediction is disabled
end
