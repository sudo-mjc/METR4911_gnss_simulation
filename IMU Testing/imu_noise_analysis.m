%% imu_noise_analysis.m
% Characterises ICM-20948 noise and bias from a STATIC log made by
% IMU_data_collection.m, and writes the results as parameters for the
% simulation (cfg.imu.* in sim_config.m).
%
% Terms estimated per axis (IEEE Std 952 Allan-variance model):
%   N  white noise          -> block NoiseDensity   (slope -1/2, read at tau = 1 s)
%   B  bias instability     -> block BiasInstability (flat floor, min / 0.664)
%   K  bias random walk     -> block RandomWalk     (slope +1/2, read at tau = 3 s)
%   b  turn-on bias         -> block ConstantBias   (mean of static data)
%
% Outputs (next to the CSV):
%   <name>_imu_meas.mat       full results struct  imu_meas
%   imu_params_measured.m     function: cfg = imu_params_measured(cfg)
%                             call it at the end of sim_config.m
%
% Needs no toolboxes (own Allan implementation). pwelch is used for a PSD
% cross-check only if the Signal Processing Toolbox is installed.

clc; clear; close all;

%% ---- User settings ----------------------------------------------------
S.file        = "";          % "" -> file picker
S.T_warmup_s  = 900;         % discard warm-up while the die heats [s]
S.lat_deg     = -27.50;      % test latitude (UQ St Lucia) for local gravity
S.h_m         = 30;          % test height above sea level [m]
S.tau_min_s   = 0.2;         % ignore tau below this (DLPF shapes short tau)
S.planar_accel_axes = [1 2]; % axes the 2D sim actually uses
S.planar_gyro_axis  = 3;

% Simulink IMU-block convention correction (from check_imu_densities, see
% eskf_build_plan.md "The sqrt(2) finding"): the block's default NoiseType
% produces N_eff = 0.707*dialog and K_eff = 1.414*dialog. To make the block
% produce the MEASURED N and K, the dialog values are compensated:
S.blk_white_comp = sqrt(2);  % dialog NoiseDensity = N_meas * sqrt(2)
S.blk_rw_comp    = 1/sqrt(2);% dialog RandomWalk   = K_meas / sqrt(2)
S.accel_sign     = -1;       % cfg.imu.accelInv: block reading = -specific force

% Datasheet values for comparison (ICM-20948 DS, one-sided densities)
DS.accel_N = 230e-6 * 9.80665;     % 230 ug/rtHz   -> m/s^2/rtHz
DS.gyro_N  = deg2rad(0.015);       % 0.015 dps/rtHz -> rad/s/rtHz

%% ---- Load -------------------------------------------------------------
if S.file == ""
    [f, p] = uigetfile("*.csv", "Select static IMU log");
    if isequal(f, 0), return; end
    S.file = string(fullfile(p, f));
end
meta = read_meta(S.file);
fs   = meta.fs_hz;  dt = 1/fs;

M = readmatrix(S.file, "CommentStyle", "#", "FileType", "text");
M = M(all(isfinite(M), 2), :);
k = M(:,1);

% ---- integrity: counter must be monotonic; fill rare gaps ----
if any(diff(k) <= 0)
    error("Sample counter is not increasing (MCU reset mid-log?). " + ...
          "Split the file at the reset and analyse one segment.");
end
n_lost = sum(diff(k) - 1);
fprintf("Samples: %d   lost: %d (%.4f %%)\n", numel(k), n_lost, 100*n_lost/(numel(k)+n_lost));
if n_lost/numel(k) > 1e-3
    warning("More than 0.1%% of samples lost - Allan results at short tau are less reliable.");
end
k_full = (k(1):k(end)).';
raw    = interp1(k, M(:,3:8), k_full, "linear");   % fills gaps only

% ---- sample-rate sanity check from the MCU tick ----
fs_tick = (numel(k)-1) / ((M(end,2) - M(1,2))/1000);
fprintf("Rate: nominal %.2f Hz, from MCU tick %.2f Hz\n", fs, fs_tick);
if abs(fs_tick - fs)/fs > 0.02
    warning("Measured rate differs from nominal by >2%% - check ODR settings.");
end

% ---- convert to SI ----
g0   = 9.80665;
acc  = raw(:,1:3) / meta.accel_lsb_per_g * g0;           % m/s^2
gyr  = deg2rad(raw(:,4:6) / meta.gyro_lsb_per_dps);      % rad/s
t    = (0:size(raw,1)-1).' * dt;

% ---- drop warm-up ----
keep = t >= S.T_warmup_s;
if nnz(keep) < 60*fs
    warning("Less than 1 min after warm-up - using all data.");
    keep(:) = true;
end
acc = acc(keep,:);  gyr = gyr(keep,:);  t = t(keep) - t(find(keep,1));
T   = t(end);
fprintf("Analysed duration: %.1f min (after %.0f s warm-up)\n\n", T/60, S.T_warmup_s);

%% ---- Static bias -------------------------------------------------------
g_loc = local_gravity(S.lat_deg, S.h_m);
mu_a  = mean(acc);
[~, gax] = max(abs(mu_a));                     % axis carrying gravity
g_vec = zeros(1,3);  g_vec(gax) = sign(mu_a(gax)) * g_loc;
bias_a = mu_a - g_vec;                         % specific-force bias, sensor frame
bias_g = mean(gyr);

tilt_deg = rad2deg(asin(min(1, norm(mu_a([1:gax-1 gax+1:3])) / g_loc)));
fprintf("Gravity on axis %d (local g = %.4f m/s^2)\n", gax, g_loc);
fprintf("Apparent tilt (if all x/y offset were tilt): %.2f deg\n", tilt_deg);
fprintf("  -> horizontal accel bias below includes tilt*g; a 6-position test separates them.\n\n");

%% ---- Allan deviation ---------------------------------------------------
L    = size(acc,1);
m    = unique(round(logspace(0, log10(floor((L-1)/2)), 150))).';
tau  = m * dt;
fprintf("Computing Allan deviation (%d cluster sizes)...\n", numel(m));
adev_a = allan_dev([acc(:,1)-mean(acc(:,1)), acc(:,2)-mean(acc(:,2)), acc(:,3)-mean(acc(:,3))], m, dt);
adev_g = allan_dev(gyr - mean(gyr), m, dt);

names = ["x","y","z"];
for i = 1:3
    A(i) = allan_fit(tau, adev_a(:,i), T, S.tau_min_s); %#ok<SAGROW>
    G(i) = allan_fit(tau, adev_g(:,i), T, S.tau_min_s); %#ok<SAGROW>
end

%% ---- PSD cross-check of N (optional) -----------------------------------
have_psd = exist("pwelch", "file") == 2;
if have_psd
    nfft = 2^nextpow2(min(L, 60*fs));
    [Pa, f] = pwelch(acc - mean(acc), hann(nfft), nfft/2, nfft, fs);
    Pg      = pwelch(gyr - mean(gyr), hann(nfft), nfft/2, nfft, fs);
    band    = f > 0.5 & f < 10;                   % flat region below the DLPF
    % pwelch is one-sided (G). White noise with Allan coefficient N has
    % two-sided PSD N^2, so N = sqrt(G/2). Datasheet densities are sqrt(G).
    N_psd_a = sqrt(median(Pa(band,:)) / 2);
    N_psd_g = sqrt(median(Pg(band,:)) / 2);
end

%% ---- Report --------------------------------------------------------------
fprintf("\n================= RESULTS (per axis x / y / z) =================\n");
fprintf("Accelerometer\n");
pr("  N  white noise    [m/s^2/rtHz]", [A.N], [A.N_ok]);
if have_psd, pr("     (PSD check)", N_psd_a, true(1,3)); end
pr("  B  bias instab.   [m/s^2]     ", [A.B], [A.B_ok]);
pr("  K  bias rand.walk [m/s^2*rtHz]", [A.K], [A.K_ok]);
pr("  b  turn-on bias   [m/s^2]     ", bias_a, true(1,3));
fprintf("  datasheet-convention density sqrt(2)*N = %s ug/rtHz  (DS: 230)\n", ...
        num2str(sqrt(2)*[A.N]/g0*1e6, "%7.1f "));
fprintf("Gyroscope\n");
pr("  N  white noise    [rad/s/rtHz]", [G.N], [G.N_ok]);
if have_psd, pr("     (PSD check)", N_psd_g, true(1,3)); end
pr("  B  bias instab.   [rad/s]     ", [G.B], [G.B_ok]);
pr("  K  bias rand.walk [rad/s*rtHz]", [G.K], [G.K_ok]);
pr("  b  turn-on bias   [rad/s]     ", bias_g, true(1,3));
fprintf("  datasheet-convention density sqrt(2)*N = %s dps/rtHz  (DS: 0.015)\n", ...
        num2str(sqrt(2)*rad2deg([G.N]), "%7.4f "));
fprintf("  ARW = %s deg/rt(h)   bias instability = %s deg/h\n", ...
        num2str(rad2deg([G.N])*60, "%6.3f "), num2str(rad2deg([G.B])*3600, "%6.2f "));
fprintf("'*' = not resolved by this record length: value is an UPPER BOUND.\n");
fprintf("Record length T = %.0f s; Allan fits use tau <= T/10 = %.0f s.\n", T, T/10);

%% ---- Package for the simulation ---------------------------------------
ia = S.planar_accel_axes;  ig = S.planar_gyro_axis;
imu_meas.source   = S.file;
imu_meas.date     = string(datetime("now"));
imu_meas.T_s      = T;
imu_meas.fs_hz    = fs;
imu_meas.settings = S;
imu_meas.accel    = struct("N", [A.N], "B", [A.B], "K", [A.K], "bias", bias_a, ...
                           "N_ok", [A.N_ok], "B_ok", [A.B_ok], "K_ok", [A.K_ok], ...
                           "tau", tau, "adev", adev_a);
imu_meas.gyro     = struct("N", [G.N], "B", [G.B], "K", [G.K], "bias", bias_g, ...
                           "N_ok", [G.N_ok], "B_ok", [G.B_ok], "K_ok", [G.K_ok], ...
                           "tau", tau, "adev", adev_g);

% Physical (effective) values for the planar sim: worst of accel x/y, gyro z.
% Taking the max keeps the filter's Q conservative on both horizontal axes.
phys.accel_N = max([A(ia).N]);   phys.accel_K = max([A(ia).K]);   phys.accel_B = max([A(ia).B]);
phys.gyro_N  = G(ig).N;          phys.gyro_K  = G(ig).K;          phys.gyro_B  = G(ig).B;
imu_meas.phys = phys;

% Block dialog values (compensated for the block's noise-type scaling)
blk.accel.NoiseDensity = phys.accel_N * S.blk_white_comp;
blk.accel.RandomWalk   = phys.accel_K * S.blk_rw_comp;
blk.accel.ConstantBias = S.accel_sign * bias_a;      % block reading = -f
blk.gyro.NoiseDensity  = phys.gyro_N  * S.blk_white_comp;
blk.gyro.RandomWalk    = phys.gyro_K  * S.blk_rw_comp;
blk.gyro.ConstantBias  = bias_g;
imu_meas.blk = blk;

[out_dir, base] = fileparts(S.file);
save(fullfile(out_dir, base + "_imu_meas.mat"), "imu_meas");
write_param_function(fullfile(out_dir, "imu_params_measured.m"), imu_meas, A, G, ia, ig);
fprintf("\nWrote %s\n  and %s\n", fullfile(out_dir, base + "_imu_meas.mat"), ...
        fullfile(out_dir, "imu_params_measured.m"));

%% ---- Plots -------------------------------------------------------------
figure("Name", "Allan deviation", "Position", [100 100 1100 450]);
plot_allan(subplot(1,2,1), tau, adev_a, A, "Accelerometer", "\sigma_A [m/s^2]", DS.accel_N/sqrt(2));
plot_allan(subplot(1,2,2), tau, adev_g, G, "Gyroscope",     "\sigma_A [rad/s]",  DS.gyro_N/sqrt(2));

figure("Name", "Static record");
win = max(1, round(10*fs));                   % 10 s moving mean shows drift
subplot(2,1,1); plot(t/60, movmean(acc - mean(acc), win)); grid on;
ylabel("Accel - mean [m/s^2]"); title("10 s moving mean (bias drift)"); legend(names);
subplot(2,1,2); plot(t/60, rad2deg(movmean(gyr - mean(gyr), win))); grid on;
ylabel("Gyro - mean [deg/s]"); xlabel("Time [min]"); legend(names);

%% ======================================================================
%  Local functions
%  ======================================================================
function meta = read_meta(file)
% Pull "# key=value" lines from the CSV header.
meta = struct("fs_hz", 100, "accel_lsb_per_g", 16384, "gyro_lsb_per_dps", 131);
fid = fopen(file, "r");  c = onCleanup(@() fclose(fid));
for i = 1:30
    s = fgetl(fid);
    if ~ischar(s) || ~startsWith(s, "#"), break; end
    tok = regexp(s, "^#\s*(\w+)=([-\d.eE+]+)\s*$", "tokens", "once");
    if ~isempty(tok) && isfield(meta, tok{1})
        meta.(tok{1}) = str2double(tok{2});
    end
end
end

function g = local_gravity(lat_deg, h)
% WGS-84 normal gravity (Somigliana) with free-air height correction.
s2 = sind(lat_deg)^2;
g  = 9.7803253359 * (1 + 0.00193185265241*s2) / sqrt(1 - 0.00669437999013*s2);
g  = g - 3.086e-6 * h;
end

function adev = allan_dev(x, m, dt)
% Overlapping Allan deviation (same estimator as MATLAB's allanvar).
% x: L-by-n rate samples, m: cluster sizes (samples), dt: sample period.
L     = size(x,1);
theta = [zeros(1,size(x,2)); cumsum(x,1)] * dt;    % integrated signal
adev  = zeros(numel(m), size(x,2));
for j = 1:numel(m)
    mj  = m(j);  tau = mj*dt;
    d   = theta(1+2*mj:L+1,:) - 2*theta(1+mj:L+1-mj,:) + theta(1:L+1-2*mj,:);
    adev(j,:) = sqrt(sum(d.^2,1) / (2*tau^2*(L-2*mj+1)));
end
end

function r = allan_fit(tau, adev, T, tau_min)
% IEEE-952 tangent-line method. Each noise term alone gives a straight line
% on the log-log plot; the measured curve can never lie BELOW any one term,
% so the tightest tangent line is an estimate that is never optimistic.
%   white:  sigma = N/sqrt(tau)          ->  N = min(sigma*sqrt(tau))
%   RW:     sigma = K*sqrt(tau/3)        ->  K = min(sigma*sqrt(3/tau))
%   flicker floor: sigma_min = 0.664*B   ->  B = min(sigma)/0.664
% A term counts as RESOLVED only if the curve actually has the matching slope
% where the tangent touches; otherwise the value is an upper bound.
v     = tau >= tau_min & tau <= T/10;     % >= ~10 clusters per point
tv    = tau(v);  av = adev(v);
slope = gradient(log10(av), log10(tv));
nv    = numel(tv);

[r.N, iN] = min(av .* sqrt(tv));
r.N_ok    = abs(slope(iN) + 0.5) < 0.15;
r.tau_N   = tv(iN);

[r.K, iK] = min(av .* sqrt(3 ./ tv));
r.K_ok    = slope(iK) > 0.35 && iK < nv - 1;
r.tau_K   = tv(iK);

[amin, iB] = min(av);
r.B       = amin / sqrt(2*log(2)/pi);   % 0.664
r.B_ok    = iB < nv - 2;                 % floor reached before the record ran out
r.tau_B   = tv(iB);
end

function s = upper_flag(ok)
if ok, s = ''; else, s = ' (UPPER BOUND - record too short)'; end
end

function pr(label, v, ok)
s = strings(1,3);
for i = 1:3
    s(i) = sprintf("%11.3e%s", v(i), repmat('*', 1, ~ok(i)));
end
fprintf("%s %s\n", label, join(s, " "));
end

function plot_allan(ax, tau, adev, R, ttl, ylab, N_ds)
axes(ax); loglog(tau, adev, "LineWidth", 1.2); hold on; grid on;
co = get(ax, "ColorOrder");
for i = 1:3
    loglog(tau, R(i).N ./ sqrt(tau), ":", "Color", co(i,:));
    loglog(tau, R(i).K .* sqrt(tau/3), "--", "Color", co(i,:));
    loglog(R(i).tau_B, R(i).B*0.664, "o", "Color", co(i,:));
end
loglog(tau, N_ds ./ sqrt(tau), "k-.");
legend(["x","y","z"], "Location", "southwest");
title(ttl + "  ( : white,  -- RW,  o floor,  -. datasheet)");
xlabel("\tau [s]"); ylabel(ylab);
end

function write_param_function(path, im, A, G, ia, ig)
% Generates a small function so sim_config.m stays the single config source:
%   cfg = imu_params_measured(cfg);   % at the end of sim_config.m
fid = fopen(path, "w");  c = onCleanup(@() fclose(fid));
w = @(varargin) fprintf(fid, varargin{:});
flag = @(ok) upper_flag(ok);
b = im.blk;
w("function cfg = imu_params_measured(cfg)\n");
w("%% IMU parameters from bench characterisation of the ICM-20948.\n");
w("%% AUTO-GENERATED by imu_noise_analysis.m on %s\n", im.date);
w("%% Source: %s\n", im.source);
w("%% Static record %.1f min at %g Hz after warm-up.\n", im.T_s/60, im.fs_hz);
w("%%\n%% Dialog values are compensated for the IMU block NoiseType scaling\n");
w("%% (white x%.4f, random walk x%.4f). Re-run check_imu_densities after\n", ...
  im.settings.blk_white_comp, im.settings.blk_rw_comp);
w("%% applying: the measured effective values should equal imu_meas.phys.\n%%\n");
w("%% Physical (effective) values used:\n");
w("%%   accel N = %.4e m/s^2/rtHz  (max of axes %s)%s\n", im.phys.accel_N, mat2str(ia), flag(all([A(ia).N_ok])));
w("%%   accel K = %.4e m/s^2*rtHz%s\n", im.phys.accel_K, flag(all([A(ia).K_ok])));
w("%%   gyro  N = %.4e rad/s/rtHz  (axis %d)%s\n", im.phys.gyro_N, ig, flag(G(ig).N_ok));
w("%%   gyro  K = %.4e rad/s*rtHz%s\n\n", im.phys.gyro_K, flag(G(ig).K_ok));

w("cfg.imu.accel.NoiseDensity = %.6e;   %% m/s^2/rtHz (dialog)\n", b.accel.NoiseDensity);
w("cfg.imu.accel.RandomWalk   = %.6e;   %% m/s^2*rtHz (dialog)\n", b.accel.RandomWalk);
w("cfg.imu.gyro.NoiseDensity  = %.6e;   %% rad/s/rtHz (dialog)\n", b.gyro.NoiseDensity);
w("cfg.imu.gyro.RandomWalk    = %.6e;   %% rad/s*rtHz (dialog)\n\n", b.gyro.RandomWalk);

w("%% Bias instability: measured, but the verified baseline keeps it at 0\n");
w("%% (a random-walk bias state cannot model flicker). Enable for the\n");
w("%% model-mismatch experiment, after checking its block scaling.\n");
w("%% cfg.imu.accel.BiasInstability = %.6e;   %% m/s^2\n", im.phys.accel_B);
w("%% cfg.imu.gyro.BiasInstability  = %.6e;   %% rad/s\n\n", im.phys.gyro_B);

w("%% Turn-on bias of THIS unit, uncalibrated, in block convention\n");
w("%% (accel sign flipped: block reading = -specific force). Use these only\n");
w("%% if the sim models NO turn-on calibration; otherwise keep the residual\n");
w("%% values. Either way, P0 bias terms must cover whatever is used here.\n");
w("%% cfg.imu.accel.ConstantBias = %s;   %% m/s^2\n", mat2str(b.accel.ConstantBias, 5));
w("%% cfg.imu.gyro.ConstantBias  = %s;   %% rad/s\n", mat2str(b.gyro.ConstantBias, 5));
w("end\n");
end
