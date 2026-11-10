function [S, R] = testt_seeds(cfg, Kp_list, Ki_list, seeds, varargin)
%TESTT_SEEDS  Closed-loop metrics for a set of (Kp, Ki) gains, averaged over noise seeds.
%
%   [S, R] = TESTT_SEEDS(cfg, Kp_list, Ki_list, seeds) runs the model once for
%   every gain pair and every seed, then reports the seed-averaged metrics.
%
%   Why several seeds: one run gives only ~40 steady-state coasts, so a few
%   coasts over sigma_thr move over_thr by several percent. Averaging over
%   seeds separates real differences between gain settings from luck.
%
%   Options (name-value):
%     'Paired'  - false (default): full grid, every Kp with every Ki.
%                 true: Kp_list(i) with Ki_list(i) only (lists must be the same length).
%     'Model'   - Simulink model name (default 'fusion_sim').
%     'Verbose' - print progress per run (default true).
%
%   Outputs:
%     S - table, one row per gain pair (seed means, plus seed std where useful).
%     R - raw metrics [nPairs x nSeeds x 5]:
%         T_mean, T_std (within run), over_thr, time_over, duty.
%
%   Example:
%     S = testt_seeds(cfg, [0.2 0.4 0.6], [0.5 0.7 0.9], 1:5);
%     heatmap(S, 'Ki', 'Kp', 'ColorVariable', 'time_over_pct');

p = inputParser;
addParameter(p, 'Paired',  false);
addParameter(p, 'Model',   'fusion_sim');
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
opt = p.Results;

% ---- Build the list of gain pairs ----------------------------------------
if opt.Paired
    assert(numel(Kp_list) == numel(Ki_list), ...
        'Paired mode needs Kp_list and Ki_list of equal length.');
    gains = [Kp_list(:), Ki_list(:)];
else
    [KP, KI] = ndgrid(Kp_list, Ki_list);
    gains = [KP(:), KI(:)];
end
nG = size(gains, 1);
nS = numel(seeds);

% Put the caller's cfg back in base when done (also on error or Ctrl-C),
% so later runs don't silently inherit the last gains/seed tested here.
cfg0    = cfg;
cleanup = onCleanup(@() restore_base(cfg0)); %#ok<NASGU>

% ---- Run ------------------------------------------------------------------
cfg.ctrl.open_loop = false;
R  = nan(nG, nS, 5);
t0 = tic;
for g = 1:nG
    cfg.ctrl.Kp = gains(g, 1);
    cfg.ctrl.Ki = gains(g, 2);
    for s = 1:nS
        cfg.sim.seed = seeds(s);
        cfg.imu.Seed = seeds(s);   cfg.gnss.Seed = seeds(s);   % set at construction, so set them here too
        assignin('base', 'cfg', cfg);                          % block dialogs read cfg
        push_params(cfg);
        out = sim(opt.Model, 'StopTime', num2str(cfg.ctrl.T_run));
        R(g, s, :) = run_metrics(out, cfg.ctrl.sigma_thr);
        if opt.Verbose
            fprintf('  [%2d/%2d] Kp=%.2f  Ki=%.2f  seed=%d   (%.0f s elapsed)\n', ...
                (g-1)*nS + s, nG*nS, gains(g,1), gains(g,2), seeds(s), toc(t0));
        end
    end
end

% ---- Summarise over seeds ---------------------------------------------------
m  = reshape(mean(R, 2, 'omitnan'),   nG, 5);
sd = reshape(std (R, 0, 2, 'omitnan'), nG, 5);

S = table(gains(:,1), gains(:,2), ...
          m(:,1), m(:,2), ...
          100*m(:,3), 100*sd(:,3), ...
          100*m(:,4), 100*sd(:,4), ...
          100*m(:,5), ...
    'VariableNames', {'Kp', 'Ki', ...
                      'T_steady', 'T_jitter', ...
                      'over_thr_pct', 'over_thr_sd', ...
                      'time_over_pct', 'time_over_sd', ...
                      'duty_pct'});

fprintf('\nClosed-loop metrics, %d seeds per setting\n', nS);
fprintf('  Kp    Ki    T_steady (s)    over_thr (coasts)   time over sigma_thr   duty\n');
for g = 1:nG
    fprintf('%4.2f  %4.2f   %5.1f +/- %3.1f    %5.1f +/- %4.1f %%     %5.2f +/- %4.2f %%     %5.2f %%\n', ...
        S.Kp(g), S.Ki(g), S.T_steady(g), S.T_jitter(g), ...
        S.over_thr_pct(g), S.over_thr_sd(g), ...
        S.time_over_pct(g), S.time_over_sd(g), S.duty_pct(g));
end
end

% ============================================================================
function M = run_metrics(out, thr)
%RUN_METRICS  [T_mean, T_std, over_thr, time_over, duty] for one run.
ce   = logical(out.coast_end.Data(:));
Tuse = out.T_sleep_active.Data(ce);   Tuse = Tuse(:);
pk   = out.sig_peak.Data(ce);         pk   = pk(:);

M = nan(1, 5);
if ~isempty(Tuse)                                   % no completed coasts -> leave NaN
    ss   = ceil(numel(Tuse)/2):numel(Tuse);         % steady half, by coast count
    M(1) = mean(Tuse(ss));
    M(2) = std(Tuse(ss));
    M(3) = mean(pk(ss) > thr);
end
sp   = out.sigma_p.Data(:);   tsp = out.sigma_p.Time;
st   = tsp >= tsp(end)/2;                           % steady half, by time
M(4) = mean(sp(st) > thr);
M(5) = mean(out.gnss_on.Data);
end

function restore_base(cfg0)
assignin('base', 'cfg', cfg0);
try
    push_params(cfg0);
catch err
    warning('testt_seeds:restore', 'Could not restore params: %s', err.message);
end
end