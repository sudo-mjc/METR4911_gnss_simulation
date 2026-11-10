%% IMU_data_collection.m
% Logs ICM-20948 data streamed by main.c over UART to a CSV file.
%
% Firmware line format (raw counts):  k,t_ms,ax,ay,az,gx,gy,gz
%   k    - sample counter (one per sensor data-ready, 100 Hz)
%   t_ms - HAL_GetTick() at read time (only used to sanity-check the rate)
%
% The CSV stores the raw counts unchanged (lossless) plus '#' comment lines
% holding the metadata needed to convert them. imu_noise_analysis.m reads it.
%
% Stop early with Ctrl+C - the file is closed cleanly and stays usable.
%
% TEST PROCEDURE (static noise run)
%   - Board on a solid, level surface, away from fans/footsteps. Don't touch it.
%   - Power on, then run this script; the analysis discards the warm-up.
%   - Duration: >= 3 h for random walk / bias instability; overnight is better.
%     A few minutes only resolves the white-noise term.

clc; clear;

%% ---- User settings ----------------------------------------------------
cfg.port       = "COM3";
cfg.baud       = 115200;
cfg.duration_s = 3*3600;          % logging time [s]; Inf = until Ctrl+C
cfg.label      = "static_flat";   % test name, goes in the filename
cfg.out_dir    = fullfile(pwd, "imu_logs");
cfg.notes      = "Bench static test, sensor z up";

% Must match the register settings in main.c
cfg.fs_hz            = 100;       % GYRO_SMPLRT_DIV = 10 -> 1100/11
cfg.accel_lsb_per_g  = 16384;     % +-2 g
cfg.gyro_lsb_per_dps = 131;       % +-250 dps
cfg.dlpf_3db_hz      = 23.9;      % DLPFCFG = 4 (accel and gyro)

%% ---- Run --------------------------------------------------------------
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
cfg.file = fullfile(cfg.out_dir, "imu_" + cfg.label + "_" + stamp + ".csv");

stats = log_imu(cfg);

fprintf("\nSaved %d samples (%.1f min) to\n  %s\n", stats.n, stats.n/cfg.fs_hz/60, cfg.file);
fprintf("Lost samples (gaps in k): %d (%.4f %%)\n", stats.lost, 100*stats.lost/max(stats.n+stats.lost,1));

%% ---- Quick look ---------------------------------------------------------
M = readmatrix(cfg.file, "CommentStyle", "#", "FileType", "text");
t = (M(:,1) - M(1,1)) / cfg.fs_hz;
acc  = M(:,3:5) / cfg.accel_lsb_per_g * 9.80665;      % m/s^2
gyro = M(:,6:8) / cfg.gyro_lsb_per_dps;               % deg/s

figure("Name", "IMU log: " + cfg.label);
subplot(2,1,1); plot(t, acc);  grid on; ylabel("Accel [m/s^2]"); legend("x","y","z");
subplot(2,1,2); plot(t, gyro); grid on; ylabel("Gyro [deg/s]");  legend("x","y","z");
xlabel("Time [s] (from sample counter)");

%% ======================================================================
function stats = log_imu(cfg)
% Reads the port in blocks (not one readline per sample) so MATLAB keeps up
% easily at 100 Hz, and writes each valid line straight to disk.

port = serialport(cfg.port, cfg.baud);
fid  = fopen(cfg.file, "w");
% onCleanup runs on normal exit AND on Ctrl+C, because this is a function.
cleaner = onCleanup(@() close_all(fid, port)); %#ok<NASGU>

% ---- metadata header (read back by imu_noise_analysis.m) ----
fprintf(fid, "# ICM-20948 raw log\n");
fprintf(fid, "# date=%s\n", string(datetime("now")));
fprintf(fid, "# label=%s\n", cfg.label);
fprintf(fid, "# notes=%s\n", cfg.notes);
fprintf(fid, "# fs_hz=%g\n", cfg.fs_hz);
fprintf(fid, "# accel_lsb_per_g=%g\n", cfg.accel_lsb_per_g);
fprintf(fid, "# gyro_lsb_per_dps=%g\n", cfg.gyro_lsb_per_dps);
fprintf(fid, "# dlpf_3db_hz=%g\n", cfg.dlpf_3db_hz);
fprintf(fid, "# columns=k,t_ms,ax,ay,az,gx,gy,gz (raw counts)\n");

flush(port);
carry  = "";
n      = 0;  lost = 0;  k_prev = NaN;
t_start = tic;  t_print = 0;

fprintf("Logging to %s\n(Ctrl+C to stop early)\n", cfg.file);
while toc(t_start) < cfg.duration_s
    nb = port.NumBytesAvailable;
    if nb == 0
        pause(0.02);
        continue
    end
    chunk = carry + string(read(port, nb, "char"));
    lines = split(chunk, newline);
    carry = lines(end);                       % last piece may be incomplete
    lines = strip(lines(1:end-1));            % strip removes the \r

    for i = 1:numel(lines)
        L = lines(i);
        v = sscanf(L, "%f,%f,%f,%f,%f,%f,%f,%f");
        if numel(v) ~= 8                      % banner / partial / corrupt line
            if strlength(L) > 0
                fprintf(fid, "# fw: %s\n", L);
                fprintf("  [fw] %s\n", L);
            end
            continue
        end
        if ~isnan(k_prev) && v(1) > k_prev + 1
            lost = lost + (v(1) - k_prev - 1);
        elseif ~isnan(k_prev) && v(1) < k_prev  % MCU reset mid-log
            fprintf(fid, "# RESET detected at row %d\n", n);
            warning("Sample counter went backwards - MCU reset?");
        end
        k_prev = v(1);
        fprintf(fid, "%s\n", L);
        n = n + 1;
    end

    if toc(t_start) - t_print > 10            % progress every 10 s
        t_print = toc(t_start);
        fprintf("  %7.0f s | %8d samples | %6.1f Hz | lost %d\n", ...
                t_print, n, n/t_print, lost);
    end
end
stats.n = n; stats.lost = lost;
end

function close_all(fid, port)
if fid > 0, fclose(fid); end
delete(port);
fprintf("Port closed, file saved.\n");
end
