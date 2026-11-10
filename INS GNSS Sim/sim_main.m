    %SIM_MAIN.M
%
% Main simulation script for data generation and processing

clear; clc;
close all;

cfg   = sim_config();
truth = gen_truth_trajectory(cfg);

t = truth.t;
q = truth.orientation;
truth_ts.position        = timeseries(truth.position,           t);
truth_ts.velocity        = timeseries(truth.velocity,           t);
truth_ts.acceleration    = timeseries(truth.acceleration,       t);
truth_ts.orientation     = timeseries(compact(q),               t);
truth_ts.angularVelocity = timeseries(truth.angularVelocity,    t);
truth_ts.psi             = timeseries(truth.psi,                t);
truth_ts.accelBody       = timeseries(truth.accelBody,          t);
%% Generate Plots
Tmax      = 360;    % s   - plotting window
v_thresh  = 0.05;   % m/s - below this counts as stationary
min_dur   = 1.0;    % s   - ignore stationary dips shorter than this

%% Crop to window
k   = truth.t <= truth.t(1) + Tmax;
t   = truth.t(k);
pos = truth.position(k,:);
psi = truth.psi(k);
acc = truth.acceleration(k,:);
spd = vecnorm(truth.velocity(k,1:2), 2, 2);

%% Detect stationary segments
isStat = spd < v_thresh;
d      = diff([0; isStat(:); 0]);
iStart = find(d == 1);
iEnd   = find(d == -1) - 1;

keep   = (t(iEnd) - t(iStart)) >= min_dur;   % drop stride-level dips
iStart = iStart(keep);
iEnd   = iEnd(keep);
nStat  = numel(iStart);
fprintf('%d stationary segments found in first %d s\n', nStat, Tmax);

%% Plot
figure('Name', 'Truth Trajectory');

% --- Path
subplot(2,2,1); hold on;
plot(pos(:,2), pos(:,1), 'LineWidth', 1.2);
for j = 1:nStat
    m = round((iStart(j) + iEnd(j))/2);              % mid-point of segment
    plot(pos(m,2), pos(m,1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k');
    text(pos(m,2), pos(m,1), sprintf('  S%d', j), 'FontSize', 8);
end
plot(pos(1,2),   pos(1,1),   'g^', 'MarkerFaceColor', 'g');   % start
plot(pos(end,2), pos(end,1), 'ks', 'MarkerFaceColor', 'k');   % end
axis equal; grid on;
xlabel('East (m)'); ylabel('North (m)');
title('a) Ground-truth path (o = stationary)');
hold off;

% --- Heading
subplot(2,2,2); hold on;
plot(t, rad2deg(psi), 'LineWidth', 1.2);
grid on; xlim([t(1) t(end)]);
xlabel('Time (s)'); ylabel('Heading (deg)');
title('b) Heading vs time');
shadeStationary(gca, t, iStart, iEnd);
hold off;

% --- Acceleration
subplot(2,2,3); hold on;
plot(t, acc(:,1), 'LineWidth', 1.2);
plot(t, acc(:,2), 'LineWidth', 1.2);
grid on; xlim([t(1) t(end)]);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('c) Acceleration (N/E components)');
legend('North', 'East', 'Location', 'southeast');
shadeStationary(gca, t, iStart, iEnd);
hold off;

% --- Speed
subplot(2,2,4); hold on;
plot(t, spd, 'LineWidth', 1.2);
grid on; xlim([t(1) t(end)]);
xlabel('Time (s)'); ylabel('Speed (m/s)');
title('d) Speed vs time');
shadeStationary(gca, t, iStart, iEnd);
hold off;

%% Simulink Simulation of INS + GNSS Integration
% Build parameter structs for running simulation blocks

prm = build_sim_params(cfg);

simInput = Simulink.SimulationInput('fusion_sim');
simInput = simInput.setVariable('cfg', cfg);

fn = fieldnames(prm);
for i = 1:numel(fn)
    simInput = simInput.setVariable(fn{i}, prm.(fn{i}));    % for the function blocks
end

% Provide timeseries data types
simInput = simInput.setVariable('truth_acc_ts',    truth_ts.acceleration);
simInput = simInput.setVariable('truth_angvel_ts', truth_ts.angularVelocity);
simInput = simInput.setVariable('truth_quat_ts',   truth_ts.orientation);
simInput = simInput.setVariable('truth_pos_ts',    truth_ts.position);
simInput = simInput.setVariable('truth_vel_ts',    truth_ts.velocity);

simInput = simInput.setModelParameter('StopTime',  num2str(cfg.sim.T), ...
                          'FixedStep', num2str(cfg.sim.dt_imu));

%% Run Simulation
out_fixed = sim(simInput);

% Pull data out of sim
%   x = [pE pN vE vN psi b_ax b_ay b_gz]'
insOLPos = out_fixed.x_ol.Data(1:2, :, :);
insOLPos = squeeze(insOLPos)';

insCLPos = out_fixed.x_nav.Data(1:2, :, :);
insCLPos = squeeze(insCLPos)';

T_plot = 90;    % Time to plot over
n = round(T_plot / cfg.sim.dt_imu);
n = min([n, size(insOLPos,1), size(insCLPos,1), size(truth.position,1)]);  % guard against short logs

figure('Name', 'INS Position Vs True Position Over time');
hold on;
plot(insOLPos(1:n,1), insOLPos(1:n,2), '-r');
plot(insCLPos(1:n,1), insCLPos(1:n,2), '-b', 'Linewidth', 1.0);
plot(truth.position(1:n,1), truth.position(1:n,2), 'k-', 'LineWidth', 0.7);
grid on;
xlabel('East [m]'); 
ylabel('North [m]');
legend('Open Loop INS ', 'Closed Loop INS', 'True Path', 'Location', 'best');
title('Truth vs. simulated INS (noise + bias) Position');
hold off;

figure('Name', 'Positional Uncertainty Over Time');
hold on;
plot(out_fixed.tout, out_fixed.sigma_p.Data);
title('Fixed Duty-Cycle Uncertainty');
yline(cfg.ctrl.sigma_set, '--', '\sigma_{set}');  yline(cfg.ctrl.sigma_thr, '-', '\sigma_{thr}');
xlabel('Time (t) [s]');
ylabel('Positional Uncertainty (\sigma_p) [m]')
grid on;
hold off;

%% Monte Carlo Simulations
run_monte = false;
if run_monte
    res = monte_carlo(truth, cfg, 30);
end

%% Incremental Controller Test
disp("Closed-loop Incremental Control")
cfg.ctrl.open_loop = false;
push_params(cfg);

out_inc  = sim('fusion_sim', 'StopTime', num2str(cfg.ctrl.T_run));
ce   = logical(out_inc.coast_end.Data(:));
Tuse = out_inc.T_sleep_active.Data(ce);      % sleep that produced each coast
cr   = logical(out_inc.crossed.Data(ce));
pk   = out_inc.sig_peak.Data(ce);
N    = numel(Tuse);   ss = ceil(N/2):N;  % steady state = second half

M.coasts   = N;
M.T_steady = [mean(Tuse(ss)) std(Tuse(ss))];
M.k_conv   = find(abs(Tuse - M.T_steady(1)) < 0.1*M.T_steady(1), 1);
M.crossed  = mean(cr(ss));               % fraction of coasts that reached sigma_set
M.over_thr = mean(pk(ss) > cfg.ctrl.sigma_thr);   % the requirement metric
M.duty     = mean(out_inc.gnss_on.Data);
disp(M)

figure('Name', 'Bang-Bang Control');  tiledlayout(2,1);
nexttile;  plot(out_inc.sigma_p.Time, out_inc.sigma_p.Data);  hold on;
title('Bang-Bang Controller Uncertainty');
yline(cfg.ctrl.sigma_set, '--', '\sigma_{set}');  yline(cfg.ctrl.sigma_thr, '-', '\sigma_{thr}');
ylabel('\sigma_p [m]');  xlabel('t [s]');  grid on;
nexttile;  plot(Tuse, 'k-', 'HandleVisibility', 'off');  hold on;
plot(find(cr),  Tuse(cr),  'o', 'DisplayName', 'Crossed');
plot(find(~cr), Tuse(~cr), '^', 'DisplayName', '~Crossed');
ylabel('T_{sleep} [s]');  xlabel('coast');  legend;  grid on;


%% PI Seeds
cfg.ctrl.open_loop = false;
cfg.ctrl.adapt_en = true;
run_seeds = false;
if run_seeds 
    disp("Closed-loop Kp/Ki Metrics");
    S = testt_seeds(cfg, [0.1 0.2 0.3], [0.2 0.3 0.4 0.5], 1:5);
    figure('Name', 'Kp/Ki sweep');
    tiledlayout(1, 2);
    nexttile;  heatmap(S, 'Ki', 'Kp', 'ColorVariable', 'time_over_pct');  title('Time over \sigma_{thr} [%]');
    nexttile;  heatmap(S, 'Ki', 'Kp', 'ColorVariable', 'over_thr_pct');   title('Coasts over \sigma_{thr} [%]');
end
%% Closed-loop run: metrics and plots.
disp("Closed-loop PI Control")
% Enable closed loop control with sleep prediction

push_params(cfg);

out_adapt  = sim('fusion_sim', 'StopTime', num2str(cfg.ctrl.T_run));
ce   = logical(out_adapt.coast_end.Data(:));
Tuse = out_adapt.T_sleep_active.Data(ce);      % sleep that produced each coast
cr   = logical(out_adapt.crossed.Data(ce));
pk   = out_adapt.sig_peak.Data(ce);
N    = numel(Tuse);   ss = ceil(N/2):N;  % steady state = second half

M.coasts   = N;
M.T_steady = [mean(Tuse(ss)) std(Tuse(ss))];
M.k_conv   = find(abs(Tuse - M.T_steady(1)) < 0.1*M.T_steady(1), 1);
M.crossed  = mean(cr(ss));               % fraction of coasts that reached sigma_set
M.over_thr = mean(pk(ss) > cfg.ctrl.sigma_thr);   % the requirement metric
M.duty     = mean(out_adapt.gnss_on.Data);
disp(M)

figure('Name', 'PI Control w/ Sleep Prediction');  tiledlayout(2,1);
nexttile;  plot(out_adapt.sigma_p.Time, out_adapt.sigma_p.Data);  hold on;
title('PI Adaptive Controller Uncertainty');
yline(cfg.ctrl.sigma_set, '--', '\sigma_{set}');  yline(cfg.ctrl.sigma_thr, '-', '\sigma_{thr}');
ylabel('\sigma_p [m]');  xlabel('t [s]');  grid on;
nexttile;  plot(Tuse, 'k-', 'HandleVisibility', 'off');  hold on;
plot(find(cr),  Tuse(cr),  'o', 'DisplayName', 'Crossed');
plot(find(~cr), Tuse(~cr), '^', 'DisplayName', '~Crossed');
ylabel('T_{sleep} [s]');  xlabel('coast');  legend;  grid on;

%% Comparison Plot: fixed baseline vs closed-loop controllers
runs   = {out_fixed, out_adapt, out_inc};
names  = ["Fixed duty-cycle", "PI", "Incremental"];
cols   = [0.25 0.25 0.25; lines(2)];      % baseline in dark grey
lstyle = ["--", "-", "-"];                % baseline dashed

figure('Name', 'Closed-loop Controller Uncertainty Comparison');
hold on;  grid on;
h = gobjects(1, numel(runs));
for k = 1:numel(runs)
    o  = runs{k};
    t  = o.sigma_p.Time;   sp = o.sigma_p.Data(:);

    % Full trace, faded so the peaks stand out
    plot(t, sp, '-', 'Color', [cols(k,:) 0.2], 'LineWidth', 0.5, ...
         'HandleVisibility', 'off');

    % Peak sigma_p at the end of each coast
    ce = logical(o.coast_end.Data(:));
    tc = o.coast_end.Time(ce);
    pk = o.sig_peak.Data(ce);   pk = pk(:);
    duty = 100 * mean(o.gnss_on.Data);
    h(k) = plot(tc, pk, lstyle(k), 'Marker', 'o', 'Color', cols(k,:), ...
                'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', cols(k,:), ...
                'DisplayName', sprintf('%s: max %.2f m, duty %.1f %%', names(k), max(pk), duty));

    % Mark the worst peak of the run
    [pmax, i] = max(pk);
    plot(tc(i), pmax, 'p', 'Color', cols(k,:), 'MarkerSize', 14, ...
         'MarkerFaceColor', cols(k,:), 'HandleVisibility', 'off');
end
yline(cfg.ctrl.sigma_set, '--', '\sigma_{set}', 'HandleVisibility', 'off');
yline(cfg.ctrl.sigma_thr, '-',  '\sigma_{thr}', 'HandleVisibility', 'off');
xline(60, ':', 'steady-state', 'HandleVisibility', 'off');
xlabel('t [s]');  ylabel('\sigma_p [m]');
title('Coast-peak uncertainty: fixed duty-cycle vs closed-loop control');
legend(h, 'Location', 'best');
hold off;
%% Power Calc
disp("Fixed Duty-Cycle");
power_fix_ctrl = gnss_power(out_fixed.gnss_on);
disp("Incremental Controller")
power_ia_ctrl = gnss_power(out_inc.gnss_on);
disp("Adaptive Controller")
power_ic_ctrl = gnss_power(out_adapt.gnss_on);

res_ctrl = monte_carlo_ctrl(truth, cfg, 30);

%% Bang Bang Controller RMSE Evaluation
disp("Bang-Bang Controller")
t   = out_inc.x_nav.Time;
X   = squeeze(out_inc.x_nav.Data);
gnss_on = out_inc.gnss_on.Data;
if size(X,1) ~= numel(t), X = X.'; end    % force N x n_states

idx_p = [1 2];            % <-- set to where pE, pN sit in x_nav
p_est = X(:, idx_p);

p_true = interp1(truth.t, truth.position(:,1:2), t);
m = eval_rmse(t, p_est, p_true, struct('t_skip', 10));          % Full run
% m = eval_rmse(t, p_est, p_true, struct('mask', ~gnss_on));      % coasting only

%% Adaptive Controller RMSE Evaluation
disp("Adaptive Controller")
t   = out_adapt.x_nav.Time;
X   = squeeze(out_adapt.x_nav.Data);
gnss_on = out_adapt.gnss_on.Data;
if size(X,1) ~= numel(t), X = X.'; end    % force N x n_states

idx_p = [1 2];            % <-- set to where pE, pN sit in x_nav
p_est = X(:, idx_p);

p_true = interp1(truth.t, truth.position(:,1:2), t);
m = eval_rmse(t, p_est, p_true, struct('t_skip', 10));          % Full run
% m = eval_rmse(t, p_est, p_true, struct('mask', ~gnss_on));      % coasting only