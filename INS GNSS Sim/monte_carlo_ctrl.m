function res = monte_carlo_ctrl(truth, cfg, N)
%MONTE_CARLO_CTRL  N-run Monte Carlo comparison of the GNSS duty-cycle strategies.
%
%   res = monte_carlo_ctrl(truth, cfg, N)        % N = 30
%
%   Runs Fixed / Incremental (bang-bang) / PI on the same truth path.
%   Run j uses IMU seed 1000+j and GNSS seed 2000+j for ALL three modes
%   (common random numbers), so each seed is a paired comparison and the
%   PI - Incremental differences are tested per seed.
%
%   Per run, per mode:
%     accuracy : rmse, cep50, cep95, emax [m], pct_in (|e| <= R_req),
%                pct_rmse_in (1 s rolling RMSE <= R_req)       [%]
%     energy   : E [J], P_avg [mW], Q [mAh], duty [%], wakes  (gnss_power)
%     control  : coasts, T_ss, T_ss_std, k_conv, crossed, over_thr
%                (same definitions as sim_main.m, steady state = 2nd half)
%
%   Logged signals needed: x_nav, gnss_on, coast_end, T_sleep_active,
%                          crossed, sig_peak
%
%   Call from sim_main.m after truth/cfg are built:
%       res_ctrl = monte_carlo_ctrl(truth, cfg, 30);

    mdl    = 'fusion_sim';
    t_skip = 10;                            % (s) as eval_rmse in sim_main
    thr    = cfg.ctrl.R_req;                % (m) accuracy threshold

    modes = struct( ...
        'name',      {'Fixed', 'Incremental', 'PI'}, ...
        'open_loop', {true,    false,         false}, ...
        'adapt_en',  {false,   false,         true});
    nM = numel(modes);

%% 1. Build N x nM simulation inputs
    t = truth.t;
    ts.acc = timeseries(truth.acceleration,       t);
    ts.ang = timeseries(truth.angularVelocity,    t);
    ts.q   = timeseries(compact(truth.orientation), t);
    ts.pos = timeseries(truth.position,           t);
    ts.vel = timeseries(truth.velocity,           t);

    n = 0;
    for j = 1:N
        for m = 1:nM
            c = cfg;
            c.imu.Seed         = 1000 + j;
            c.gnss.Seed        = 2000 + j;
            c.ctrl.open_loop   = modes(m).open_loop;
            c.ctrl.adapt_en    = modes(m).adapt_en;
            p = build_sim_params(c);

            n = n + 1;
            in(n) = Simulink.SimulationInput(mdl);          %#ok<AGROW>
            in(n) = in(n).setVariable('cfg', c);
            fn = fieldnames(p);
            for i = 1:numel(fn), in(n) = in(n).setVariable(fn{i}, p.(fn{i})); end
            in(n) = in(n).setVariable('truth_acc_ts',    ts.acc);
            in(n) = in(n).setVariable('truth_angvel_ts', ts.ang);
            in(n) = in(n).setVariable('truth_quat_ts',   ts.q);
            in(n) = in(n).setVariable('truth_pos_ts',    ts.pos);
            in(n) = in(n).setVariable('truth_vel_ts',    ts.vel);
            in(n) = in(n).setModelParameter('StopTime',  num2str(cfg.sim.T), ...
                                            'FixedStep', num2str(cfg.sim.dt_imu));
        end
    end
    out = parsim(in, 'ShowProgress', 'on');   % no PCT: for n=1:numel(in), out(n)=sim(in(n)); end

%% 2. Metrics per run
    n = 0;
    for j = 1:N
        for m = 1:nM
            n = n + 1;
            o = out(n);
            if ~isempty(o.ErrorMessage)
                error('Run %d (%s, seed %d) failed:\n%s', n, modes(m).name, j, o.ErrorMessage);
            end
            r = struct();

            % --- Accuracy
            [X, tx] = getsig(o, 'x_nav');
            pT = interp1(truth.t, truth.position(:,1:2), tx);
            k  = tx >= t_skip;
            e  = X(k,1:2) - pT(k,:);
            e2 = sum(e.^2, 2);
            en = sqrt(e2);
            win = max(1, round(1.0 / cfg.sim.dt_imu));           % 1 s window
            r.rmse        = sqrt(mean(e2));
            r.cep50       = quant(en, 0.50);
            r.cep95       = quant(en, 0.95);
            r.emax        = max(en);
            r.pct_in      = 100 * mean(en <= thr);
            r.pct_rmse_in = 100 * mean(sqrt(movmean(e2, win)) <= thr);

            % --- Energy (your gnss_power, printout suppressed)
            gp = [];
            evalc('gp = gnss_power(o.gnss_on);');
            r.E     = pick(gp, {'E','E_J','energy'});
            r.P_avg = pick(gp, {'P_avg','Pavg','P_mW'});
            r.Q     = pick(gp, {'Q','Q_mAh','charge'});
            r.wakes = pick(gp, {'nWake','wakes','n_wakes','nWakes'});
            r.duty  = 100 * mean(getsig(o, 'gnss_on'));

            % --- Controller (definitions from sim_main.m)
            ce   = getsig(o, 'coast_end') > 0.5;
            Tsl  = getsig(o, 'T_sleep_active');   Tuse = Tsl(ce);
            cr   = getsig(o, 'crossed') > 0.5;    cr   = cr(ce);
            pk   = getsig(o, 'sig_peak');         pk   = pk(ce);
            nc   = numel(Tuse);   ss = ceil(nc/2):nc;
            r.coasts   = nc;
            r.T_ss     = mean(Tuse(ss));
            r.T_ss_std = std(Tuse(ss));
            kc = find(abs(Tuse - r.T_ss) < 0.1*r.T_ss, 1);
            if isempty(kc), kc = NaN; end
            r.k_conv   = kc;
            r.crossed  = mean(cr(ss));
            r.over_thr = mean(pk(ss) > cfg.ctrl.sigma_thr);

            R(j, m) = r;                                     %#ok<AGROW>
        end
    end

%% 3. Statistics
    flds  = {'rmse','cep50','cep95','emax','pct_in','pct_rmse_in', ...
             'E','P_avg','Q','duty','wakes', ...
             'coasts','T_ss','T_ss_std','k_conv','crossed','over_thr'};
    stats = struct();
    for f = flds
        V = reshape([R.(f{1})], N, nM);                      % N x nM
        stats.(f{1}).data = V;
        for m = 1:nM
            [mu, sd, ci] = msci(V(:,m));
            stats.(f{1}).mean(m) = mu;
            stats.(f{1}).std(m)  = sd;
            stats.(f{1}).ci(m,:) = ci;
        end
        [mu, ~, ci] = msci(V(:,3) - V(:,2));                 % paired PI - Incremental
        stats.(f{1}).d_pi_inc    = mu;
        stats.(f{1}).d_pi_inc_ci = ci;
    end

    E = stats.E.data;
    for m = 2:nM
        s    = 1 - E(:,m) ./ E(:,1);                         % saving per seed
        life = 1 ./ (1 - s) - 1;                             % battery-life increase
        [res.save_mean(m), res.save_std(m)] = msci(100*s);
        [res.life_mean(m), res.life_std(m)] = msci(100*life);
    end

    res.modes = {modes.name};
    res.N     = N;
    res.runs  = R;
    res.stats = stats;
    res.pass  = sum(stats.pct_rmse_in.data >= 100*cfg.ctrl.p_req, 1);

%% 4. Report
    L = repmat('-', 1, 100);
    fprintf('\n%s\n Duty-cycle controller Monte Carlo: %d runs per mode, mean +/- std [95%% CI of mean]\n%s\n', L, N, L);
    fprintf('%-13s', 'metric');  fprintf('%-29s', res.modes{:});  fprintf('PI - Inc (paired)\n');
    for f = flds
        S = stats.(f{1});
        fprintf('%-13s', f{1});
        for m = 1:nM
            fprintf('%-29s', sprintf('%.4g +/- %.3g [%.4g, %.4g]', ...
                    S.mean(m), S.std(m), S.ci(m,1), S.ci(m,2)));
        end
        sig = S.d_pi_inc_ci(1) > 0 || S.d_pi_inc_ci(2) < 0;
        fprintf('%+.3g [%+.3g, %+.3g] %s\n', S.d_pi_inc, S.d_pi_inc_ci, ...
                ifelse(sig, '*', ''));
    end
    fprintf('%s\n (* = paired 95%% CI excludes zero)\n\n', L);
    for m = 2:nM
        fprintf(' %-12s vs Fixed: energy saving %5.1f +/- %.1f %%   GNSS battery life +%5.1f +/- %.1f %%\n', ...
                res.modes{m}, res.save_mean(m), res.save_std(m), res.life_mean(m), res.life_std(m));
    end
    fprintf('\n');
    for m = 1:nM
        fprintf(' %-12s goal (>= %.0f %% within %.2f m): passed %d/%d runs, worst %.1f %%\n', ...
                res.modes{m}, 100*cfg.ctrl.p_req, thr, res.pass(m), N, ...
                min(stats.pct_rmse_in.data(:,m)));
    end
    fprintf('%s\n\n', L);

%% 5. Plot: energy and RMSE per run, coloured by mode
    show = {'P_avg', 'Average GNSS power [mW]'; ...
            'rmse', 'Horizontal RMSE [m]'};
    cols = [0.85 0.33 0.10;      % Fixed       - orange
            0.00 0.45 0.74;      % Incremental - blue
            0.47 0.67 0.19];     % PI          - green
    jit  = 0.12 * (rand(N,1) - 0.5);                         % same jitter in both tiles

    figure('Name', 'MC duty-cycle controllers');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    for i = 1:size(show,1)
        nexttile;  hold on;  grid on;  box on;
        V = stats.(show{i,1}).data;                          % N x nM
        S = stats.(show{i,1});
        for m = 1:nM
            % individual runs
            scatter(m + jit, V(:,m), 22, cols(m,:), 'filled', ...
                    'MarkerFaceAlpha', 0.45, 'HandleVisibility', 'off');
            % mean +/- 95% CI
            errorbar(m + 0.28, S.mean(m), S.mean(m) - S.ci(m,1), S.ci(m,2) - S.mean(m), ...
                     'o', 'Color', cols(m,:) * 0.7, 'MarkerFaceColor', cols(m,:), ...
                     'MarkerSize', 8, 'LineWidth', 1.8, 'CapSize', 10, ...
                     'DisplayName', sprintf('%s: %.3g \\pm %.2g', res.modes{m}, S.mean(m), S.std(m)));
        end
        xlim([0.5 nM+0.6]);  xticks(1:nM);  xticklabels(res.modes);
        ylabel(show{i,2});
        legend('Location', 'best');
    end
    sgtitle(sprintf('%d-run Monte Carlo: individual runs (dots) and mean \\pm 95%% CI', N));
end

% =====================================================================
function [mu, sd, ci] = msci(x)
%MSCI  Mean, std and 95% CI of the mean (NaNs ignored, no toolboxes).
    x  = x(~isnan(x));
    mu = mean(x);  sd = std(x);
    if numel(x) < 2, ci = [NaN NaN]; return; end
    hw = tcrit(numel(x) - 1) * sd / sqrt(numel(x));
    ci = [mu - hw, mu + hw];
end

function t = tcrit(df)
    tab = [12.71 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 ...
           2.201 2.179 2.160 2.145 2.131 2.120 2.110 2.101 2.093 2.086 ...
           2.080 2.074 2.069 2.064 2.060 2.056 2.052 2.048 2.045 2.042];
    if df <= numel(tab), t = tab(df); else, t = 1.96; end
end

function q = quant(x, p)
    x = sort(x(:));
    q = x(max(1, ceil(p * numel(x))));
end

function v = pick(s, names)
%PICK  First matching field of the gnss_power output.
    for i = 1:numel(names)
        if isstruct(s) && isfield(s, names{i}), v = double(s.(names{i})); return; end
    end
    error(['monte_carlo_ctrl: gnss_power output has no field named %s.\n' ...
           'Fields returned: %s\nAdd the right name to the pick() call.'], ...
           strjoin(names, '/'), strjoin(fieldnames(s), ', '));
end

function y = ifelse(c, a, b)
    if c, y = a; else, y = b; end
end

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