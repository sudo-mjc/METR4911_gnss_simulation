function res = monte_carlo_v2(truth, cfg, N)
%MONTE_CARLO_V2  N-run Monte Carlo NEES / NIS consistency test of the ESKF.
%
%   res = monte_carlo_v2(truth, cfg, N)        % N = 30
%
%   Each run uses new IMU and GNSS noise seeds on the same truth path.
%     NEES_k = e_k' P_k^-1 e_k   e = position error, P = position block   (2 dof)
%     NIS_k  = y_k' S_k^-1 y_k   y = full [pos; vel] innovation, from the
%                                eskf_update_v2 'nis' output, fix ticks (4 dof)
%
%   The N-run average of a chi-square(k) variable, times N, is chi-square(kN),
%   so the 95% acceptance band for the averages is
%
%     [chi2inv(0.025, kN), chi2inv(0.975, kN)] / N
%       NEES (k = 2, N = 30): [1.35, 2.78]     NIS (k = 4, N = 30): [3.05, 5.07]
%
%   The filter passes if ~95% of the averaged values lie inside the band.
%   Run with a FIXED duty cycle (cfg.ctrl.open_loop = true) so every run
%   has its GNSS fixes on the same ticks and the NIS can be averaged.
%
%   Logged signals needed: x_nav, P_diag, P_xy, nis, gnss_valid

    t_settle = 60;                          % (s) ignore the initial convergence
    mdl      = 'fusion_sim';

%% 1. Run N simulations
    for j = 1:N
        c = cfg;
        c.imu.Seed  = 1000 + j;
        c.gnss.Seed = 2000 + j;
        p = build_sim_params(c);

        in(j) = Simulink.SimulationInput(mdl);              %#ok<AGROW>
        in(j) = in(j).setVariable('cfg', c);
        fn = fieldnames(p);
        for i = 1:numel(fn), in(j) = in(j).setVariable(fn{i}, p.(fn{i})); end
        in(j) = in(j).setVariable('truth_acc_ts',    truth.acceleration);
        in(j) = in(j).setVariable('truth_angvel_ts', truth.angularVelocity);
        in(j) = in(j).setVariable('truth_quat_ts',   truth.orientation);
        in(j) = in(j).setVariable('truth_pos_ts',    truth.position);
        in(j) = in(j).setVariable('truth_vel_ts',    truth.velocity);
        in(j) = in(j).setModelParameter('StopTime',  num2str(cfg.sim.T), ...
                                        'FixedStep', num2str(cfg.sim.dt_imu));
    end
    out = parsim(in, 'ShowProgress', 'on');   % no PCT: for j=1:N, out(j)=sim(in(j)); end

%% 2. NEES and NIS per run, summed across runs
    for j = 1:N
        [xnav, t] = getsig(out(j), 'x_nav');
        Pd  = getsig(out(j), 'P_diag');
        Pxy = getsig(out(j), 'P_xy');
        nis = getsig(out(j), 'nis');                              % 4-dof NIS from eskf_update_v2
        fix = getsig(out(j), 'gnss_valid') > 0.5 & nis > 0;       % update ticks

        if j == 1
            pT       = interp1(truth.t, truth.position(:,1:2), t);
            k        = t >= t_settle;
            nees_sum = zeros(size(t));
            nis_sum  = zeros(size(t));
            nis_cnt  = zeros(size(t));
            e2_sum   = 0;
        end

        % NEES (2x2 inverse written out)
        e    = xnav(:,1:2) - pT;
        nees = (Pd(:,2).*e(:,1).^2 - 2*Pxy.*e(:,1).*e(:,2) + Pd(:,1).*e(:,2).^2) ...
               ./ (Pd(:,1).*Pd(:,2) - Pxy.^2);

        nis(~fix) = 0;

        nees_sum = nees_sum + nees;
        nis_sum  = nis_sum  + nis;
        nis_cnt  = nis_cnt  + fix;
        e2_sum   = e2_sum   + sum(sum(e(k,:).^2, 2));
    end

%% 3. Averages and acceptance band
    % 95% band for the N-run mean of chi2(dof): chi2inv([.025 .975], dof*N) / N
    chi2band = @(dof) 2*gammaincinv([0.025 0.975], dof*N/2) / N;
    band     = chi2band(2);                           % NEES: position, 2 dof
    band_nis = chi2band(4);                           % NIS: position + velocity, 4 dof

    nees_bar = nees_sum(k) / N;
    t_nees   = t(k);

    kf      = k & nis_cnt == N;                       % fix ticks shared by all runs
    nis_bar = nis_sum(kf) / N;
    t_nis   = t(kf);

    res.band         = band;
    res.band_nis     = band_nis;
    res.nees_mean    = mean(nees_bar);
    res.nees_inband  = 100 * mean(nees_bar >= band(1) & nees_bar <= band(2));
    res.nis_mean     = mean(nis_bar);
    res.nis_inband   = 100 * mean(nis_bar  >= band_nis(1) & nis_bar  <= band_nis(2));
    res.n_fix        = numel(nis_bar);
    res.rmse         = sqrt(e2_sum / (nnz(k) * N));
    res.t_nees = t_nees;  res.nees_bar = nees_bar;
    res.t_nis  = t_nis;   res.nis_bar  = nis_bar;

    if res.n_fix == 0
        warning('No fix ticks are shared by all runs: use a fixed duty cycle for the NIS test.');
    end

%% 4. Report
    L = repmat('-', 1, 62);
    fprintf('\n%s\n ESKF consistency, %d-run Monte Carlo, t >= %g s\n%s\n', L, N, t_settle, L);
    fprintf('                    dof  ideal   95%% band        mean   inside band\n');
    fprintf(' NEES (pos)          2    2.0   [%.2f, %.2f]   %6.3f   %6.1f %%\n', ...
            band, res.nees_mean, res.nees_inband);
    fprintf(' NIS  (pos + vel)    4    4.0   [%.2f, %.2f]   %6.3f   %6.1f %%   (%d fix epochs)\n', ...
            band_nis, res.nis_mean, res.nis_inband, res.n_fix);
    fprintf(' target: >= 95 %% of the %d-run means inside the band\n', N);
    fprintf('%s\n Position RMSE  %.3f m\n%s\n\n', L, res.rmse, L);

%% 5. Plot
    figure('Name', 'MC NEES / NIS');
    subplot(2,1,1);
    plot(t_nees, nees_bar, 'b-', 'LineWidth', 0.6); hold on;
    yline(2, 'k-', 'LineWidth', 1.2);
    yline(band(1), 'r--');  yline(band(2), 'r--');
    grid on;  ylim([0 4]);  ylabel(sprintf('NEES (%d-run mean)', N));
    title(sprintf('NEES: mean %.2f, %.1f%% inside 95%% band', res.nees_mean, res.nees_inband));

    subplot(2,1,2);
    plot(t_nis, nis_bar, 'b.', 'MarkerSize', 4); hold on;
    yline(4, 'k-', 'LineWidth', 1.2);
    yline(band_nis(1), 'r--');  yline(band_nis(2), 'r--');
    grid on;  ylim([0 8]);  xlabel('Time [s]');  ylabel(sprintf('NIS, 4 dof (%d-run mean)', N));
    title(sprintf('NIS: mean %.2f, %.1f%% inside 95%% band', res.nis_mean, res.nis_inband));
end

% =====================================================================
function [v, t] = getsig(data, name)
%GETSIG  Pull a logged signal from either To Workspace or signal logging.
    ts = [];
    try
        ts = data.(name);
    catch
    end
    if isempty(ts)
        ts = data.logsout.getElement(name).Values;
    end
    t = ts.Time(:);
    v = squeeze(ts.Data);
    if size(v,1) ~= numel(t) && size(v,2) == numel(t)
        v = v.';
    end
    v = double(v);
end