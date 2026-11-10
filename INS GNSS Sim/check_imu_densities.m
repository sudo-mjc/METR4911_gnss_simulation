function eff = check_imu_densities(cfg, truth)
%CHECK_IMU_DENSITIES  Measure the noise densities the IMU model actually produces.
%
%   eff = check_imu_densities(cfg, truth)
%
%   The value you type into NoiseDensity / RandomWalk is not necessarily the
%   density the block generates -- MathWorks' NoiseType scaling has differed
%   between releases, and the ratio is typically 1.000, 0.707 or 1.414.
%   Since the ESKF builds Qc directly from these numbers, a scaling mismatch
%   shows up as a uniformly optimistic or conservative NEES with a correct NIS.
%
%   This uses the imuSensor System object directly -- the Simulink IMU block
%   wraps the same object, so it takes seconds and needs no model run.
%
%   Each noise term is isolated by constructing a sensor with ONLY that term
%   enabled (accelparams/gyroparams default every other term to zero).
%
%   Returns the MEASURED effective densities. Feed these into eskf_prm:
%       eskf_prm.Na = eff.Na;   eskf_prm.Ng = eff.Ng;
%       eskf_prm.Ka = eff.Ka;   eskf_prm.Kg = eff.Kg;

    dt = cfg.sim.dt_imu;
    fs = 1/dt;

    Na_set = cfg.imu.accel.VelocityRandomWalk;       % (m/s^2)/sqrt(Hz)
    Ng_set = cfg.imu.gyro.AngleRandomWalk;           % (rad/s)/sqrt(Hz)
    Ka_set = cfg.imu.accel.AccelerationRandomWalk;   % (m/s^2)/sqrt(s)
    Kg_set = cfg.imu.gyro.RateRandomWalk;            % (rad/s)/sqrt(s)

    %% ---- white measurement noise only --------------------------------
    [ra, rg] = run_sensor(cfg, truth, fs, ...
        accelparams('NoiseDensity', Na_set*[1 1 1]), ...
        gyroparams( 'NoiseDensity', Ng_set*[1 1 1]));

    % white noise: sigma_discrete = N / sqrt(dt)   ->   N = sigma * sqrt(dt)
    eff.Na = mean(std(ra)) * sqrt(dt);
    eff.Ng =      std(rg)  * sqrt(dt);

    %% ---- bias random walk only ---------------------------------------
    [ba, bg] = run_sensor(cfg, truth, fs, ...
        accelparams('RandomWalk', Ka_set*[1 1 1]), ...
        gyroparams( 'RandomWalk', Kg_set*[1 1 1]));

    % with no white noise the residual IS the bias trajectory;
    % its increments have std K*sqrt(dt)  ->  K = std(diff(b)) / sqrt(dt)
    eff.Ka = mean(std(diff(ba))) / sqrt(dt);
    eff.Kg =      std(diff(bg))  / sqrt(dt);

    %% ---- report -------------------------------------------------------
    rows = { 'Na  accel white ', Na_set, eff.Na
             'Ng  gyro  white ', Ng_set, eff.Ng
             'Ka  accel bias RW', Ka_set, eff.Ka
             'Kg  gyro  bias RW', Kg_set, eff.Kg };

    fprintf('\n%s\n', repmat('-',1,62));
    fprintf(' IMU effective noise densities (measured from the block)\n');
    fprintf('%s\n', repmat('-',1,62));
    fprintf(' %-18s %12s %12s %8s\n', 'term', 'set', 'measured', 'ratio');
    for i = 1:size(rows,1)
        ratio = rows{i,3} / rows{i,2};
        fprintf(' %-18s %12.4e %12.4e %8.3f%s\n', ...
                rows{i,1}, rows{i,2}, rows{i,3}, ratio, flag(ratio));
    end
    fprintf('%s\n', repmat('-',1,62));
    fprintf([' ratio 1.000 -> the dialog value IS the density; Q is already right\n' ...
             ' ratio 0.707 -> block is quieter than set; your Q is 2x too large\n' ...
             ' ratio 1.414 -> block is louder than set; your Q is 2x too small\n' ...
             ' Either way: set eskf_prm.Na/Ng/Ka/Kg from the MEASURED column.\n\n']);
end

% ---------------------------------------------------------------------
function s = flag(r)
    if     abs(r-1.000) < 0.05, s = '   (1)';
    elseif abs(r-0.707) < 0.05, s = '   (1/sqrt2)';
    elseif abs(r-1.414) < 0.05, s = '   (sqrt2)';
    else,                       s = '   <-- unexpected';
    end
end

% ---------------------------------------------------------------------
function [res_a, res_g] = run_sensor(cfg, truth, fs, ap, gp)
%RUN_SENSOR  Residual of the sensor output against the analytic truth.
    imu = imuSensor('accel-gyro', ...
        'SampleRate',     fs, ...
        'ReferenceFrame', 'ENU', ...
        'Accelerometer',  ap, ...
        'Gyroscope',      gp, ...
        'Seed',           cfg.imu.Seed);

    [a, g] = imu(truth.acceleration, truth.angularVelocity, truth.orientation);

    % accelInv puts the reading into the filter's convention, where it should
    % equal the analytic body specific force plus whatever noise is enabled.
    res_a = cfg.imu.accelInv * a(:,1:2) - truth.accelBody(:,1:2);
    res_g = g(:,3)                      - truth.angularVelocity(:,3);
end
